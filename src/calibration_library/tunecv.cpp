#include "tunecv.h"

void readPersonalisationConductionVelocity (std::string filename, std::map<std::string,double> &cv_map) {
    char buffer[256];
    FILE *file = fopen(filename.c_str(),"r");
    if (!file) {
        printf("[-] ERROR! File '%s' not found!\n", filename.c_str());
        exit(EXIT_FAILURE);
    }
    uint32_t num_lines = 0;
    const char *delimiters = " \t\n\r";
    // Read each line of the file
    while (fgets(buffer, sizeof(buffer), file) != NULL) {
        // Build a vector with the tokens using whitespaces as delimiter
        char *token = strtok(buffer, delimiters);
        std::vector<std::string> tokens;
        while (token != NULL) {
            tokens.push_back(std::string(token));
            token = strtok(NULL, delimiters);
        }    
        // Skip line 0
        if (num_lines > 0) {
            if (tokens.size() == 4) {
                char *end;
                cv_map[tokens[0]] = strtod(tokens[2].c_str(), &end);
            }
            else {
                printf("[-] ERROR! Invalid number of tokens!\n");
                exit(EXIT_FAILURE);
            }
        }
        num_lines++;
    } 
    fclose(file);
}

double calculate_conduction_velocity_from_simulation ()
{
    char filename[500];
    sprintf(filename, "%s/outputs/tunecv/cable/tissue_activation_time_map_pulse_it_0.vtu", CONVERTER_PATH);

    // Read all the data from the file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename);
    reader->Update();

    vtkUnstructuredGrid* unstructuredGrid = reader->GetOutput();
    uint32_t num_points = unstructuredGrid->GetNumberOfPoints();
    uint32_t num_lines = unstructuredGrid->GetNumberOfCells();

    vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
    cellLocator->SetDataSet(unstructuredGrid);
    cellLocator->BuildLocator();

    // Points of interest
    double point_x_0[3] = {23000, 100, 100};
    double point_x_1[3] = {27000, 100, 100};

    // Find (closest points) Cell indexes for CV computation
    vtkIdType cellId_x_0; // the cell id of the cell containing the closest point
    vtkIdType cellId_x_1; // the cell id of the cell containing the closest point
    
    double closestPoint_x_0[3];
    double closestPoint_x_1[3];
    
    int subId; // not needed
    double closestPointDist2; // not needed
    
    cellLocator->FindClosestPoint(point_x_0, closestPoint_x_0, cellId_x_0, subId, closestPointDist2);
    cellLocator->FindClosestPoint(point_x_1, closestPoint_x_1, cellId_x_1, subId, closestPointDist2);
    
    double delta_s_x = sqrt(pow(closestPoint_x_0[0]-closestPoint_x_1[0], 2)); // only changes over x-axis
    
    cout << delta_s_x << endl;
    
    // Read points scalar values
    std::string array_name = "Scalars_";
    vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast(unstructuredGrid->GetCellData()->GetArray(array_name.c_str()));
    double cv_x = -1.0;

    if(array) {        
        double delta_lat_x = (array->GetValue(cellId_x_1) - array->GetValue(cellId_x_0)); // ms
        
	    cout << delta_lat_x << endl;
	    cout << array->GetValue(cellId_x_0) << endl;
	    cout << array->GetValue(cellId_x_1) << endl;
	
        cv_x = (delta_s_x / delta_lat_x)*0.001;     // {m/s}
	printf("%g", cv_x);
    }
    else {
        cerr << "[!] ERROR! No 'Scalar_value' found for the points!" << endl;
        exit(EXIT_FAILURE);
    }
    return cv_x;
}

void write_configuration_file (const double sigma) {
    char filename[500];
    sprintf(filename, "%s/outputs/tunecv/configs/cable.ini", CONVERTER_PATH);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"[main]\n");
    fprintf(file,"num_threads=12\n");
    fprintf(file,"dt_pde=0.02\n");
    fprintf(file,"simulation_time=2000.0\n");
    fprintf(file,"abort_on_no_activity=false\n");
    fprintf(file,"use_adaptivity=false\n");
    fprintf(file,"quiet=true\n");
    fprintf(file,"\n");
    
    fprintf(file,"[update_monodomain]\n");
    fprintf(file,"main_function=update_monodomain_default\n");
    fprintf(file,"library_file=%s/shared_libs/libdefault_update_monodomain.so\n", MONOALG_PATH);
    fprintf(file,"\n");
    
    // For saving the LATs in a format that can be read for calculating the CVs
    fprintf(file,"[save_result]\n");
    fprintf(file,"print_rate=1\n");
    fprintf(file,"output_dir=%s/outputs/tunecv/cable\n", CONVERTER_PATH);
    fprintf(file,"save_pvd=true\n");
    fprintf(file,"file_prefix=V\n");
    fprintf(file,"save_activation_time=true\n");
    fprintf(file,"activation_threshold_tissue=0.3\n");
    fprintf(file,"save_apd=false\n");
    fprintf(file,"library_file=%s/shared_libs/libdefault_save_mesh_purkinje.so\n", MONOALG_PATH);
    fprintf(file,"main_function=save_tissue_with_activation_times\n");
    fprintf(file,"init_function=init_save_tissue_with_activation_times\n");
    fprintf(file,"end_function=end_save_tissue_with_activation_times\n");
    fprintf(file,"remove_older_simulation=true\n");
    fprintf(file,"\n");
    
    fprintf(file,"[assembly_matrix]\n");
    fprintf(file,"init_function=set_initial_conditions_fvm\n");
    fprintf(file,"sigma_x=%g\n",sigma);
    fprintf(file,"sigma_y=%g\n",sigma);
    fprintf(file,"sigma_z=%g\n",sigma);
    fprintf(file,"library_file=%s/shared_libs/libdefault_matrix_assembly.so\n", MONOALG_PATH);
    fprintf(file,"main_function=homogeneous_sigma_assembly_matrix\n");
    fprintf(file,"\n");
    
    fprintf(file,"[linear_system_solver]\n");
    fprintf(file,"tolerance=1e-16\n");
    fprintf(file,"use_preconditioner=no\n");
    fprintf(file,"max_iterations=500\n");
    fprintf(file,"library_file=%s/shared_libs/libdefault_linear_system_solver.so\n", MONOALG_PATH);
    fprintf(file,"main_function=conjugate_gradient\n");
    fprintf(file,"use_gpu=no\n");
    fprintf(file,"\n");
    
    fprintf(file,"[domain]\n");  
    fprintf(file,"name=Simple Cable\n");
    fprintf(file,"start_dx=200.0\n");
    fprintf(file,"start_dy=200.0\n");
    fprintf(file,"start_dz=200.0\n");
    fprintf(file,"cable_length=50000.0\n");
    fprintf(file,"library_file=%s/shared_libs/libdefault_domains.so\n", MONOALG_PATH);
    fprintf(file,"main_function=initialize_grid_with_cable_mesh\n");
    fprintf(file,"\n");
    
    fprintf(file,"[ode_solver]\n");
    fprintf(file,"dt=0.02\n");
    fprintf(file,"use_gpu=yes\n");
    fprintf(file,"gpu_id=0\n");
    fprintf(file,"library_file=%s/shared_libs/libminimal_model.so\n", MONOALG_PATH);
    fprintf(file,"\n");
    
    fprintf(file,"[stim_benchmark]\n");
    fprintf(file,"start = 0.0\n");
    fprintf(file,"duration = 2.0\n");
    fprintf(file,"current = 1.0\n");
    fprintf(file, "min_x = 0.0\n");
    fprintf(file, "max_x = 3000.0\n");
    fprintf(file, "min_y = 0.0\n");
    fprintf(file, "max_y = 3000.0\n");
    fprintf(file, "min_z = 0.0\n");
    fprintf(file, "max_z = 3000.0\n");
    fprintf(file,"main_function=stim_x_y_z_limits\n");
    fprintf(file,"library_file=%s/shared_libs/libdefault_stimuli.so\n", MONOALG_PATH);
    fprintf(file,"\n");

    fclose(file);
}

double runTuneCV (double target_cv) {
    double cv, factor;
    double sigma = 0.0002;

    do {
        write_configuration_file(sigma);
        
        // Run the simulation
        char command[500];
        sprintf(command,"%s/bin/MonoAlg3D -c %s/outputs/tunecv/configs/cable.ini", MONOALG_PATH, CONVERTER_PATH);
        system(command);
        
        cv = calculate_conduction_velocity_from_simulation();
        factor = pow(target_cv/cv,2);
        sigma = sigma*factor;

        printf("\n|| Target CV = %g m/s || Computed CV = %g m/s || Factor = %g || Adjusted sigma = %g mS/um ||\n\n",target_cv,cv,factor,sigma);

    }while ( fabs(cv-target_cv) > TOLERANCE );
    
    printf("\n[+] Target conductivity = %.10lf mS/um\n",sigma);
    return sigma;
}

void calibrateConductionVelocity(std::map<std::string,double> cv_map, std::map<std::string,double> &conductivity_map) {
    
    for (auto it = cv_map.begin(); it != cv_map.end(); ++it) {
        if (!std::isnan(it->second)) {
            printf("==================================================================\n");
            printf("[!] Calibrating %s = %lf ...\n", it->first.c_str(), it->second);
            printf("==================================================================\n");
            double target_cv = it->second;
            double target_sigma = runTuneCV(target_cv);
            conductivity_map[it->first] = target_sigma;
        }
    }
}