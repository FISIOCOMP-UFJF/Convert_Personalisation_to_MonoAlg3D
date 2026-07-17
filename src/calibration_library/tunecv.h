#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkHexahedron.h>
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkCellLocator.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkTransformPolyDataFilter.h>

static const double TOLERANCE = 1.0e-03;
static const char MONOALG_PATH[500] = "/home/berg/Github/versions_MonoAlg3D_C/daniel_moreira/MonoAlg3D_C";
static const char CONVERTER_PATH[500] = "/home/berg/MEGA/PostDoc/CNPq/Conhecimento_Brasil_2024/Programas/Convert_Personalisation_to_MonoAlg3D";

void readPersonalisationConductionVelocity (std::string filename, std::map<std::string,double> &cv_map);
void calibrateConductionVelocity(std::map<std::string,double> cv_map, std::map<std::string,double> &conductivity_map);
void write_stimulus_monoalg3d_format (vtkSmartPointer<vtkFloatArray> lat_array, \
            vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter);
void write_configuration_file (const double sigma);
double runTuneCV (double target_cv);

#endif