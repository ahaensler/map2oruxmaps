/**********************************************************************************************
    Copyright (C) 2012 Oliver Eichler oliver.eichler@gmx.de
                  2015 Adrian Hänsler

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

**********************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <list>
#include <wctype.h>
#include <errno.h>
#include <cpl_conv.h>

#include <gdalwarper.h>
#include <gdal_priv.h>
#include <proj_api.h>
#include <ogr_spatialref.h>
#include <sqlite3.h>
#include <sys/stat.h>

#include "CInputFile.h"

#define WHAT_STR "map2oruxmaps, Version 1.0"

using namespace std;



static bool sort_in_files(CInputFile& f1, CInputFile& f2) {
    return f1.getXScale() < f2.getXScale();
}

void db_exec(sqlite3 *db, char *sql){
    char *err_msg = 0;
    int rc = sqlite3_exec(db, sql, 0, 0, &err_msg);
    if (rc != SQLITE_OK ) {
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);        
        sqlite3_close(db);
        exit(-1);
    }
}



int main(int argc, char ** argv)
{
    uint32_t nLevels  =  0;
    int quality       = -1;
    int subsampling   = -1;
    int crop_to_base  =  0;
    int skip_next_arg =  0;
    string mapname;
    string basepath;
    string dbfile;
    string xmlfile;
    string gdal_cache_setting = "128";
    GDALRIOResampleAlg resampling_alg = GRIORA_NearestNeighbour;
    list<string> infilenames;
    list<CInputFile> infiles;


    printf("\n****** %s ******\n", WHAT_STR);

    if(argc < 2) {
      fprintf(stderr,"\nusage: map2oruxmaps -r <algo> -q <1..100> -s <411|422|444> -m <memory> -c -n <mapname> <file1> <file2> ... <fileN> <output directory>\n");
      fprintf(stderr,"\n");
      fprintf(stderr,"  -r    Resampling algorithm (nearest, bilinear, cubic, cubicspline, lanczos, average, mode, gauss) \n");
      fprintf(stderr,"  -q    The JPEG quality from 1 to 100. Default is 75\n");
      fprintf(stderr,"  -s    The chroma subsampling. Default is 411\n");
      fprintf(stderr,"  -m    Sets GDAL_CACHEMAX memory (for big maps).\n");
      fprintf(stderr,"  -c    Crop all zoom levels to the base map.\n");
      fprintf(stderr,"\n");
      fprintf(stderr,"\n");
      exit(-1);
    }

    GDALAllRegister();

    for(int i = 1; i < (argc - 1); i++) {
      if (skip_next_arg) {
        skip_next_arg--;
        continue;
      }

      if (argv[i][0] == '-') {
        if (towupper(argv[i][1]) == 'N') {
          mapname.assign( argv[i+1] );
          skip_next_arg = 1;
          continue;
        } else if (towupper(argv[i][1]) == 'Q') {
          quality = atol(argv[i+1]);
          skip_next_arg = 1;
          continue;
        } else if (towupper(argv[i][1]) == 'S') {
          subsampling = atol(argv[i+1]);
          skip_next_arg = 1;
          continue;
        } else if (towupper(argv[i][1]) == 'C') {
          crop_to_base = 1;
          continue;
        } else if (towupper(argv[i][1]) == 'M') {
          gdal_cache_setting = argv[i+1];
          skip_next_arg = 1;
          continue;
        } else if (towupper(argv[i][1]) == 'R'){
          if(!strcmp(argv[i+1], "nearest"))
            resampling_alg = GRIORA_NearestNeighbour;
          else if(!strcmp(argv[i+1], "bilinear"))
            resampling_alg = GRIORA_Bilinear;
          else if(!strcmp(argv[i+1], "cubic"))
            resampling_alg = GRIORA_Cubic;
          else if(!strcmp(argv[i+1], "cubicspline"))
            resampling_alg = GRIORA_CubicSpline;
          else if(!strcmp(argv[i+1], "lanczos"))
            resampling_alg = GRIORA_Lanczos;
          else if(!strcmp(argv[i+1], "average"))
            resampling_alg = GRIORA_Average;
          else if(!strcmp(argv[i+1], "mode"))
            resampling_alg = GRIORA_Mode;
          else if(!strcmp(argv[i+1], "gauss"))
            resampling_alg = GRIORA_Gauss;
          else {
            fprintf(stderr,"Unknown resampling algorithm.\n");
            exit(-1);
          }
          skip_next_arg = 1;
          continue;
        }
      }
      infilenames.push_back( argv[i] );
    }

    CPLSetConfigOption( "GDAL_CACHEMAX", gdal_cache_setting.c_str() );

    basepath.assign( argv[argc-1] );
    if(basepath[basepath.length()-1]!='/') basepath += "/";
    if(mapname==""){
      fprintf(stderr,"A map name is required.\n");
      exit(-1);
    }
    basepath += mapname;
    if(mkdir(basepath.c_str(), 0666)){
      if(errno != EEXIST){
        fprintf(stderr,"Could not create the map directory. %s\n", basepath.c_str());
        exit(-1);
      }
    }
    for(auto fn : infilenames){
      CInputFile cif(fn, TILESIZE, mapname, resampling_alg);
      infiles.push_back( cif );
    }
    basepath += "/";
    dbfile = basepath + "OruxMapsImages.db";
    xmlfile = basepath + mapname + ".otrk2.xml";


    infiles.sort(sort_in_files);

    if(crop_to_base)
      for(auto it=next(infiles.begin()); it!=infiles.end(); ++it)
        it->crop_to(infiles.front());


    double xscale = infiles.front().getXScale();
    double yscale = infiles.front().getYScale();
    for(auto it=infiles.begin(); it!=prev(infiles.end()); ++it)
      nLevels += it->calcLevels(next(it)->getXScale(), xscale, yscale);
    nLevels += infiles.back().calcLevels(0.0, xscale, yscale);

    for(auto i : infiles)
      i.summarize();

    printf("\n\nThere is a total of %i tiles to process.", CInputFile::getTilesTotal());
    printf("\n\n");


    sqlite3 *db;
    
    unlink(dbfile.c_str());
    int rc = sqlite3_open(dbfile.c_str(), &db);
    if (rc != SQLITE_OK) {
      fprintf(stderr, "Cannot open database: %s\n", sqlite3_errmsg(db));
      sqlite3_close(db);
      return 1;
    }
    

    const char sql[] = 
                "CREATE TABLE android_metadata (locale TEXT);"
                "CREATE TABLE tiles (x int, y int, z int, image blob, PRIMARY KEY (x,y,z));"
                "CREATE INDEX IND on tiles (x,y,z);"
                "INSERT INTO android_metadata VALUES('en_US');";
    db_exec(db, (char*)sql);
    

    ofstream xmlstream;
    xmlstream.open(xmlfile.c_str(), ios::out);

    xmlstream << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      <<"<OruxTracker xmlns:orux=\"http://oruxtracker.com/app/res/calibration\"\n"
      <<" versionCode=\"3.0\">\n"
      <<"<MapCalibration layers=\"true\" layerLevel=\"0\">\n"
      <<"<MapName><![CDATA["<<mapname.c_str()<<"]]></MapName>\n";


    // write layers
    sqlite3_exec(db, "BEGIN", NULL, NULL, NULL);
    for(auto i:infiles)
      i.writeLevels(db, xmlstream, quality, subsampling);
    sqlite3_exec(db, "END", NULL, NULL, NULL);
    sqlite3_close(db);


    xmlstream<<"</MapCalibration>\n</OruxTracker>";


    GDALDestroyDriverManager();
    printf("\n");

    return 0;
}
