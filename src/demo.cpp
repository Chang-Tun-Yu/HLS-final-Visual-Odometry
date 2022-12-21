/*
Copyright 2012. All rights reserved.
Institute of Measurement and Control Systems
Karlsruhe Institute of Technology, Germany

This file is part of libviso2.
Authors: Andreas Geiger

libviso2 is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or any later version.

libviso2 is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
libviso2; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301, USA 
*/

/*
  Documented C++ sample code of stereo visual odometry (modify to your needs)
  To run this demonstration, download the Karlsruhe dataset sequence
  '2010_03_09_drive_0019' from: www.cvlibs.net!
  Usage: ./viso2 path/to/sequence/2010_03_09_drive_0019
*/

#include <iostream>
#include <string>
#include <vector>
#include <stdint.h>
#include <fstream>

#include <viso_stereo.h>
#include <viso_mono.h>
#include <png++/png.hpp>

using namespace std;

int main_mono (int argc, char** argv) {

  // we need the path name to 2010_03_09_drive_0019 as input argument
  if (argc<2) {
    cerr << "Usage: ./viso2 path/to/sequence/2010_03_09_drive_0019" << endl;
    return 1;
  }

  // sequence directory
  string dir = argv[1];
  
  // set most important visual odometry parameters
  // for a full parameter list, look at: viso_stereo.h
  // VisualOdometryStereo::parameters param;
  VisualOdometryMono::parameters param;
  
  // calibration parameters for sequence 2010_03_09_drive_0019 
  // param.calib.f  = 645.24; // focal length in pixels
  // param.calib.cu = 635.96; // principal point (u-coordinate) in pixels
  // param.calib.cv = 194.13; // principal point (v-coordinate) in pixels
  // param.base     = 0.5707; // baseline in meters
  param.calib.f  = 491.543832; // focal length in pixels
  param.calib.cu = 484.474328; // principal point (u-coordinate) in pixels
  param.calib.cv = 147.888234; // principal point (v-coordinate) in pixels
  param.height   = 1.6;
  param.pitch    = -0.08;


  
  // init visual odometry
  // VisualOdometryStereo viso(param);
  VisualOdometryMono viso(param);
  
  // current pose (this matrix transforms a point from the current
  // frame's camera coordinates to the first frame's camera coordinates)
  Matrix pose = Matrix::eye(4);
    
  // loop through all frames i=0:372
  for (int32_t i=0; i<373; i++) {

    // input file names
    // char base_name[256]; sprintf(base_name,"%06d.png",i);
    char binary_name[256]; sprintf(binary_name,"%06d.dat",i);
    // string img_file_name  = dir + "/I1_" + base_name;
    string img_binary_name  = dir + "/I1_" + binary_name;
    
    // catch image read/write errors here
    try {

      // load left and right input image
      // png::image< png::gray_pixel > Img(img_file_name);

      // image dimensions
      // int32_t width  = Img.get_width();
      // int32_t height = Img.get_height();
      int32_t width  = 1024;
      int32_t height = 284;

      // convert input images to uint8_t buffer
      uint8_t* Img_data  = (uint8_t*)malloc(width*height*sizeof(uint8_t));
      // int32_t k=0;
      // for (int32_t v=0; v<height; v++) {
      //   for (int32_t u=0; u<width; u++) {
      //     Img_data[k]  = Img.get_pixel(u,v);
      //     k++;
      //   }
      // }
      ifstream ifs;
      ifs.open(img_binary_name.c_str(), ifstream::binary | ifstream::in);
      ifs.read((char*) Img_data, width*height*sizeof(uint8_t));
      ifs.close();

      // status
      cout << "Processing: Frame: " << i;
      
      // compute visual odometry
      int32_t dims[] = {width,height,width};
      if (viso.process(Img_data,dims)) {
      
        // on success, update current pose
        pose = pose * Matrix::inv(viso.getMotion());
      
        // output some statistics
        double num_matches = viso.getNumberOfMatches();
        double num_inliers = viso.getNumberOfInliers();
        cout << ", Matches: " << num_matches;
        cout << ", Inliers: " << 100.0*num_inliers/num_matches << " %" << ", Current pose: " << endl;
        cout << pose << endl << endl;

      } else {
        cout << " ... failed!" << endl;
      }

      // release uint8_t buffers
      free(Img_data);

    // catch image read errors here
    } catch (...) {
      cerr << "ERROR: Couldn't read input files!" << endl;
      return 1;
    }
  }
  
  // output
  cout << "Demo complete! Exiting ..." << endl;

  // exit
  return 0;
}

int main_stereo (int argc, char** argv) {

  // we need the path name to 2010_03_09_drive_0019 as input argument
  if (argc<2) {
    cerr << "Usage: ./viso2 path/to/sequence/2010_03_09_drive_0019" << endl;
    return 1;
  }

  // sequence directory
  string dir = argv[1];
  
  // set most important visual odometry parameters
  // for a full parameter list, look at: viso_stereo.h
  VisualOdometryStereo::parameters param;
  
  // calibration parameters for sequence 2010_03_09_drive_0019 
  param.calib.f  = 491.543832; // focal length in pixels
  param.calib.cu = 484.474328; // principal point (u-coordinate) in pixels
  param.calib.cv = 147.888234; // principal point (v-coordinate) in pixels
  param.base     = 0.5707; // baseline in meters
  
  // init visual odometry
  VisualOdometryStereo viso(param);
  
  // current pose (this matrix transforms a point from the current
  // frame's camera coordinates to the first frame's camera coordinates)
  Matrix pose = Matrix::eye(4);
    
  // loop through all frames i=0:372
  for (int32_t i=0; i<373; i++) {

    // input file names
    // char base_name[256]; sprintf(base_name,"%06d.png",i);
    char binary_name[256]; sprintf(binary_name,"%06d.dat",i);
    // string left_img_file_name  = dir + "/I1_" + base_name;
    // string right_img_file_name = dir + "/I2_" + base_name;
    string left_binary_file_name  = dir + "/I1_" + binary_name;
    string right_binary_file_name = dir + "/I2_" + binary_name;
    
    // catch image read/write errors here
    try {

      // load left and right input image
      // png::image< png::gray_pixel > left_img(left_img_file_name);
      // png::image< png::gray_pixel > right_img(right_img_file_name);

      // image dimensions
      // int32_t width  = left_img.get_width();
      // int32_t height = left_img.get_height();
      int32_t width  = 1024;
      int32_t height = 284;

      // convert input images to uint8_t buffer
      uint8_t* left_img_data  = (uint8_t*)malloc(width*height*sizeof(uint8_t));
      uint8_t* right_img_data = (uint8_t*)malloc(width*height*sizeof(uint8_t));
      // int32_t k=0;
      // for (int32_t v=0; v<height; v++) {
      //   for (int32_t u=0; u<width; u++) {
      //     left_img_data[k]  = left_img.get_pixel(u,v);
      //     right_img_data[k] = right_img.get_pixel(u,v);
      //     k++;
      //   }
      // }
      ifstream ifs;
      ifs.open(left_binary_file_name.c_str(), ifstream::binary | ifstream::in);
      ifs.read((char*) left_img_data, width*height*sizeof(uint8_t));
      ifs.close();
      ifs.open(right_binary_file_name.c_str(), ifstream::binary | ifstream::in);
      ifs.read((char*) right_img_data, width*height*sizeof(uint8_t));
      ifs.close();
      // status
      cout << "Processing: Frame: " << i;
      
      // compute visual odometry
      int32_t dims[] = {width,height,width};
      if (viso.process(left_img_data,right_img_data,dims)) {
      
        // on success, update current pose
        pose = pose * Matrix::inv(viso.getMotion());
      
        // output some statistics
        double num_matches = viso.getNumberOfMatches();
        double num_inliers = viso.getNumberOfInliers();
        cout << ", Matches: " << num_matches;
        cout << ", Inliers: " << 100.0*num_inliers/num_matches << " %" << ", Current pose: " << endl;
        cout << pose << endl << endl;

      } else {
        cout << " ... failed!" << endl;
      }

      // saving binary img files
      // ofstream ofileleft(left_binary_file_name.c_str(), ofstream::binary | ofstream::out);
      // ofileleft.write((char*)left_img_data, width*height*sizeof(uint8_t));
      // ofileleft.close();
      // ofstream ofileright(right_binary_file_name.c_str(), ofstream::binary | ofstream::out);
      // ofileright.write((char*)right_img_data, width*height*sizeof(uint8_t));
      // ofileright.close();

      // release uint8_t buffers
      free(left_img_data);
      free(right_img_data);
      
    // catch image read errors here
    } catch (...) {
      cerr << "ERROR: Couldn't read input files!" << endl;
      return 1;
    }
    
  }
  
  // output
  cout << "Demo complete! Exiting ..." << endl;

  // exit
  return 0;
}

int main (int argc, char** argv) {
  bool mono = true;
  if (mono) {
    cout << "========== [MONO VO] ==========" << endl;
    main_mono(argc, argv);
  }
  else {
    cout << "========== [STEREO VO] ==========" << endl;
    main_stereo(argc, argv);
  }
  return 0;
}