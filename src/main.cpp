/*
  Author: @bergolho (Lucas Arantes Berg)
  Last update: 17/07/2026

  Program that converts the output from the reduced-order Eikonal Personalisation simulation to a 
  biophysically monodomain one.
  
  1) Apply a Transform filter over the root nodes with LAT from the "Cardiac_Personalisation_Brazil"
    program. The transformation consist only of a scaling from {cm} to {um}.
  2) Run the 'tuneCV' method to calibrate the inferred conduction velocities from the Personalisation to
    monodomain conductivities.
    Reference: 
    Costa, C. M., Hoetzl, E., Rocha, B. M., Prassl, A. J. & Plank, G. Automatic parameterization strategy for cardiac electrophysiology simulations. In Computing in Cardiology 2013 373–376 (IEEE, 2013).

    Inputs: 
      - Biventricular Personalisation LAT map
      - Inferred Personalisation root nodes with LAT
      - Inferred Personalisation conduction velocities
    
    Outputs:
      - Inferred Personalisation root nodes with LAT in micrometers for MonoAlg3D
      - Calibrated conductivities for MonoAlg3D
      - Stimulus section of the root nodes for MonoAlg3D
  */

#include <iostream>
#include <string>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkTransformPolyDataFilter.h>

#include "calibration_library/tunecv.h"

using namespace std;

void write_matrix_assembly_monoalg3d_format(std::map<std::string,double> conductivity_map) {
  
  char filename[500];
  sprintf(filename,"%s/outputs/monoalg3d_conductivities.txt",CONVERTER_PATH);
  FILE *file = fopen(filename,"w+");

  for (auto it = conductivity_map.begin(); it != conductivity_map.end(); ++it) {
    if (it->first == "fibre_speed") {
      fprintf(file,"sigma_l = %.10lf\n", it->second);
    }
    else if (it->first == "sheet_speed") {
      fprintf(file,"sigma_t = %.10lf\n", it->second);
    }
    else if (it->first == "normal_speed") {
      fprintf(file,"sigma_n = %.10lf\n", it->second);
    }
    else if (it->first == "endo_dense_speed") {
      fprintf(file,"sigma_l_dense = %.10lf\n", it->second);
      fprintf(file,"sigma_t_dense = %.10lf\n", it->second);
      fprintf(file,"sigma_n_dense = %.10lf\n", it->second);
    }
    else if (it->first == "endo_sparse_speed") {
      fprintf(file,"sigma_l_sparse = %.10lf\n", it->second);
      fprintf(file,"sigma_t_sparse = %.10lf\n", it->second);
      fprintf(file,"sigma_n_sparse = %.10lf\n", it->second);
    }
    else if (it->first == "purkinje_speed") {
      fprintf(file,"sigma_purkinje = %.10lf\n", it->second);
    }
    else if (it->first == "fibrosis_speed") {
      fprintf(file,"sigma_fibrosis = %.10lf\n", it->second);
    }
    else {
      printf("[-] ERROR! Invalid conduction velocity '%s'!\n", it->first.c_str());
      exit(EXIT_FAILURE);
    }
  }
}

void write_stimulus_monoalg3d_format (vtkSmartPointer<vtkFloatArray> lat_array, \
            vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter) {
  
  char filename[500];
  sprintf(filename,"%s/outputs/monoalg3d_stimulus_section.txt",CONVERTER_PATH);
  FILE *file = fopen(filename,"w+");

  const double stim_duration = 2.0;
  const double stim_current = 1.0;
  const double stim_radius = 500.0;
  
  uint32_t num_root_nodes = transformFilter->GetOutput()->GetNumberOfPoints();
  for (uint32_t i = 0; i < num_root_nodes; i++) {
    fprintf(file,"[stim_%d]\n", i+1);
    fprintf(file,"start = %.2lf\n", lat_array->GetValue(i));
    fprintf(file,"duration = %lf\n", stim_duration);
    fprintf(file,"current = %lf\n", stim_current);
    fprintf(file,"center_x = %lf\n", transformFilter->GetOutput()->GetPoint(i)[0]);
    fprintf(file,"center_y = %lf\n", transformFilter->GetOutput()->GetPoint(i)[1]);
    fprintf(file,"center_z = %lf\n", transformFilter->GetOutput()->GetPoint(i)[2]);
    fprintf(file,"radius = %lf\n", stim_radius);
    fprintf(file,"main_function = stim_sphere\n\n");
  }
  fclose(file);
}

void transformRootNodes (std::string input_biv_filename, std::string input_root_nodes_filename, std::string output_root_nodes_filename) {
  // 1) Read UnstructuredGrid biventricular mesh from Personalisation and discover the (x,y,z) bounds
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(input_biv_filename.c_str());
  reader->Update();

  vtkUnstructuredGrid *unstructured_grid = reader->GetOutput();

  double bounds[6];
  unstructured_grid->GetBounds(bounds);

  double xmin = bounds[0];
  double xmax = bounds[1];
  double ymin = bounds[2];
  double ymax = bounds[3];
  double zmin = bounds[4];
  double zmax = bounds[5];

  // 2) Read Polydata root nodes with LAT
  vtkSmartPointer<vtkGenericDataObjectReader> reader_root_nodes = vtkSmartPointer<vtkGenericDataObjectReader>::New();
  reader_root_nodes->SetFileName(input_root_nodes_filename.c_str());
  reader_root_nodes->Update();

  // All of the standard data types can be checked and obtained like this:
  vtkPolyData *polydata_grid = NULL;
  if (reader_root_nodes->IsFilePolyData()) {
    std::cout << "Output is polydata," << std::endl;
    polydata_grid = reader_root_nodes->GetPolyDataOutput();
    std::cout << "   output has " << polydata_grid->GetNumberOfPoints() << " points." << std::endl;
    
    std::string array_name = "LAT";
    vtkSmartPointer<vtkFloatArray> lat_array = vtkFloatArray::SafeDownCast(polydata_grid->GetPointData()->GetArray(array_name.c_str()));

    // 3) Apply a Transform filter.
    // Translate to the origin with the geometry bounds
    // Scale from centimeters to micrometers
    const double cm_to_um = 1e+4;
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Translate(-xmin*cm_to_um,-ymin*cm_to_um,-zmin*cm_to_um);
    transform->Scale(cm_to_um, cm_to_um, cm_to_um);

    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();

    transformFilter->SetInputData(polydata_grid);
    transformFilter->SetTransform(transform);

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(output_root_nodes_filename.c_str());
    writer->SetInputConnection(transformFilter->GetOutputPort());
    writer->Write();

    // 4) Output in MonoAlg3D format
    write_stimulus_monoalg3d_format(lat_array,transformFilter);
  }
}

int main (int argc, char *argv[]) {
  if (argc-1 != 4) {
    printf("Usage:> %s <input_biv_lat_file> <input_root_nodes_file> <input_inferred_cvs_file> <output_root_nodes_file>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  // Get input data
  std::string input_biv_filename = argv[1];
  std::string input_root_nodes_filename = argv[2];
  std::string input_inferred_cvs_filename = argv[3];
  std::string output_root_nodes_filename = argv[4];

  transformRootNodes(input_biv_filename, input_root_nodes_filename, output_root_nodes_filename);

  std::map<std::string,double> cv_map;
  readPersonalisationConductionVelocity(input_inferred_cvs_filename, cv_map);

  std::map<std::string,double> conductivity_map;
  calibrateConductionVelocity(cv_map, conductivity_map);

  write_matrix_assembly_monoalg3d_format(conductivity_map);

  return 0;
}
