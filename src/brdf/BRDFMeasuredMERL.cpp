/*
Copyright Disney Enterprises, Inc. All rights reserved.

This license governs use of the accompanying software. If you use the software, you
accept this license. If you do not accept the license, do not use the software.

1. Definitions
The terms "reproduce," "reproduction," "derivative works," and "distribution" have
the same meaning here as under U.S. copyright law. A "contribution" is the original
software, or any additions or changes to the software. A "contributor" is any person
that distributes its contribution under this license. "Licensed patents" are a
contributor's patent claims that read directly on its contribution.

2. Grant of Rights
(A) Copyright Grant- Subject to the terms of this license, including the license
conditions and limitations in section 3, each contributor grants you a non-exclusive,
worldwide, royalty-free copyright license to reproduce its contribution, prepare
derivative works of its contribution, and distribute its contribution or any derivative
works that you create.
(B) Patent Grant- Subject to the terms of this license, including the license
conditions and limitations in section 3, each contributor grants you a non-exclusive,
worldwide, royalty-free license under its licensed patents to make, have made,
use, sell, offer for sale, import, and/or otherwise dispose of its contribution in the
software or derivative works of the contribution in the software.

3. Conditions and Limitations
(A) No Trademark License- This license does not grant you rights to use any
contributors' name, logo, or trademarks.
(B) If you bring a patent claim against any contributor over patents that you claim
are infringed by the software, your patent license from such contributor to the
software ends automatically.
(C) If you distribute any portion of the software, you must retain all copyright,
patent, trademark, and attribution notices that are present in the software.
(D) If you distribute any portion of the software in source code form, you may do
so only under this license by including a complete copy of this license with your
distribution. If you distribute any portion of the software in compiled or object code
form, you may only do so under a license that complies with this license.
(E) The software is licensed "as-is." You bear the risk of using it. The contributors
give no express warranties, guarantees or conditions. You may have additional
consumer rights under your local laws which this license cannot change.
To the extent permitted under your local laws, the contributors exclude the
implied warranties of merchantability, fitness for a particular purpose and non-
infringement.
*/

#include <cstdlib>
#include <string>
#include <fstream>
#include <string.h>
#include <regex>
#include <cstdlib>
#include "BRDFMeasuredMERL.h"
#include "DGLShader.h"
#include "Paths.h"

#define BRDF_SAMPLING_RES_THETA_H       90
#define BRDF_SAMPLING_RES_THETA_D       90
#define BRDF_SAMPLING_RES_PHI_D         360



BRDFMeasuredMERL::BRDFMeasuredMERL()
                 : brdfData(nullptr), paramsPack(nullptr),
				   cmbBRDFData(nullptr), cmbParamsPack(nullptr), cmbTBO(NULL), cmbEnabled(false)
{
    std::string path = (getShaderTemplatesPath() + "measured.func");

    // read the shader
    std::ifstream ifs( path.c_str() );
    std::string temp;
    while( getline( ifs, temp ) )
        brdfFunction += (temp + "\n");
}



BRDFMeasuredMERL::~BRDFMeasuredMERL()
{
	delete[] brdfData;
	if (paramsPack != nullptr) {
		delete paramsPack;
	}
	if (cmbBRDFData != nullptr) {
		delete[] cmbBRDFData;
	}
	if (cmbParamsPack != nullptr) {
		delete cmbParamsPack;
	}

    glBindBuffer(GL_TEXTURE_BUFFER_EXT, tbo);
    glDeleteBuffers( 1, &tbo);

	if (cmbTBO != NULL) {
		glBindBuffer(GL_TEXTURE_BUFFER_EXT, cmbTBO);
		glDeleteBuffers(1, &cmbTBO);
	}
}


std::string BRDFMeasuredMERL::getBRDFFunction()
{
    return brdfFunction;
}


bool BRDFMeasuredMERL::loadMERLData( const char* filename )
{
    // the BRDF's name is just the filename
    name = std::string(filename);

    // read in the MERL BRDF data
    FILE *f = fopen(filename, "rb");
    if (!f)
            return false;
    
    int dims[3];
    if (fread(dims, sizeof(int), 3, f) != 3) {
        fprintf(stderr, "read error\n");
        fclose(f);
        return false;
    }
    numBRDFSamples = dims[0] * dims[1] * dims[2];
    if (numBRDFSamples != BRDF_SAMPLING_RES_THETA_H *
                            BRDF_SAMPLING_RES_THETA_D *
                            BRDF_SAMPLING_RES_PHI_D / 2) 
    {
        fprintf(stderr, "Dimensions don't match\n");
        fclose(f);
        return false;
    }
    
    // read the data
    double* brdf = (double*) malloc (sizeof(double)*3*numBRDFSamples);
    if (fread(brdf, sizeof(double), 3*numBRDFSamples, f) != size_t(3*numBRDFSamples)) {
        fprintf(stderr, "read error\n");
        fclose(f);
	return false;
    }
    fclose(f);

    // now transform it to RGBA floats
    brdfData = new float[ numBRDFSamples * 3 ];
    for( int i = 0; i < numBRDFSamples; i++ )
    {
            brdfData[i*3 + 0] = brdf[i*3 + 0];
            brdfData[i*3 + 1] = brdf[i*3 + 1];
            brdfData[i*3 + 2] = brdf[i*3 + 2];
    }

    // now we can dump the old data
    free( brdf );

	// read pre-computed parameters
	std::fstream paramsfs("mat_params.txt", std::ios::in); // TODO: hard-coded file name
	if (!paramsfs.is_open()) {
		fprintf(stderr, "can not read pre-computed parameters file\n");
	}
	else {
		auto p1 = name.find_last_of('/') + 1;
		auto p2 = name.find_last_of('.');
		auto mat_name = name.substr(p1, p2 - p1);
		std::string params;
		std::regex re(mat_name +
					  R"(\s([0-9]*\.?[0-9]+)\s\((\S+),(\S+),(\S+)\)\s\((\S+),(\S+),(\S+)\)\s\((\S+),(\S+),(\S+)\)\s(\d+)(\**))");
		while (std::getline(paramsfs, params)) {
			std::smatch matches;
			if (std::regex_match(params, matches, re) && paramsPack == nullptr) {
				// TODO: potential duplicate storage
				paramsPack = new MERLParametersPack();
				paramsPack->ggx = std::atof(matches[1].str().c_str());
				paramsPack->diffuseAlbedo[0] = std::atof(matches[2].str().c_str());
				paramsPack->diffuseAlbedo[1] = std::atof(matches[3].str().c_str());
				paramsPack->diffuseAlbedo[2] = std::atof(matches[4].str().c_str());
				paramsPack->specularAlbedo[0] = std::atof(matches[5].str().c_str());
				paramsPack->specularAlbedo[1] = std::atof(matches[6].str().c_str());
				paramsPack->specularAlbedo[2] = std::atof(matches[7].str().c_str());
				paramsPack->optimalThreshold[0] = std::atof(matches[8].str().c_str());
				paramsPack->optimalThreshold[1] = std::atof(matches[9].str().c_str());
				paramsPack->optimalThreshold[2] = std::atof(matches[10].str().c_str());
				paramsPack->cluster = std::atoi(matches[11].str().c_str());

				addFloatParameter("diffuseAlbedo_R", 0, 1, std::atof(matches[2].str().c_str()));
				addFloatParameter("diffuseAlbedo_G", 0, 1, std::atof(matches[3].str().c_str()));
				addFloatParameter("diffuseAlbedo_B", 0, 1, std::atof(matches[4].str().c_str()));
				addFloatParameter("specularAlbedo_R", 0, 1, std::atof(matches[5].str().c_str()));
				addFloatParameter("specularAlbedo_G", 0, 1, std::atof(matches[6].str().c_str()));
				addFloatParameter("specularAlbedo_B", 0, 1, std::atof(matches[7].str().c_str()));
			}
		}
	}

    return true;
}

bool BRDFMeasuredMERL::loadCmbMERL(const char* filename) {
	// read in the MERL BRDF data
	FILE *f = fopen(filename, "rb");
	if (!f)
		return false;

	int dims[3];
	if (fread(dims, sizeof(int), 3, f) != 3) {
		fprintf(stderr, "read error\n");
		fclose(f);
		return false;
	}
	numBRDFSamples = dims[0] * dims[1] * dims[2];
	if (numBRDFSamples != BRDF_SAMPLING_RES_THETA_H *
		BRDF_SAMPLING_RES_THETA_D *
		BRDF_SAMPLING_RES_PHI_D / 2) {
		fprintf(stderr, "Dimensions don't match\n");
		fclose(f);
		return false;
	}

	// read the data
	double* brdf = (double*)malloc(sizeof(double) * 3 * numBRDFSamples);
	if (fread(brdf, sizeof(double), 3 * numBRDFSamples, f) != size_t(3 * numBRDFSamples)) {
		fprintf(stderr, "read error\n");
		fclose(f);
		return false;
	}
	fclose(f);

	bool firstLoad = false;
	// now transform it to RGBA floats
	if (cmbBRDFData == nullptr) {
		firstLoad = true;
		cmbBRDFData = new float[numBRDFSamples * 3];
	}
	for (int i = 0; i < numBRDFSamples; i++) {
		cmbBRDFData[i * 3 + 0] = brdf[i * 3 + 0];
		cmbBRDFData[i * 3 + 1] = brdf[i * 3 + 1];
		cmbBRDFData[i * 3 + 2] = brdf[i * 3 + 2];
	}

	// now we can dump the old data
	free(brdf);

	unsigned int numBytes = numBRDFSamples * 3 * sizeof(float);
	if (firstLoad) {
		glGenBuffers(1, &cmbTBO);
		glBindBuffer(GL_TEXTURE_BUFFER_EXT, cmbTBO);
		glBufferData(GL_TEXTURE_BUFFER_EXT, numBytes, 0, GL_STATIC_DRAW);

		glGenTextures(1, &cmbTEX);
		glBindTexture(GL_TEXTURE_BUFFER_EXT, cmbTEX);
		glTexBufferEXT(GL_TEXTURE_BUFFER_EXT, GL_INTENSITY32F_ARB, cmbTBO);
		glBindBuffer(GL_TEXTURE_BUFFER_EXT, 0);
	}
	glBindBuffer(GL_TEXTURE_BUFFER_EXT, cmbTBO);
	float* p = (float*)glMapBuffer(GL_TEXTURE_BUFFER_EXT, GL_WRITE_ONLY);
	memcpy(p, cmbBRDFData, numBytes);
	glUnmapBuffer(GL_TEXTURE_BUFFER_EXT);
	glBindBuffer(GL_TEXTURE_BUFFER_EXT, 0);

	// read pre-computed parameters
	std::fstream paramsfs("mat_params.txt", std::ios::in);
	if (!paramsfs.is_open()) {
		fprintf(stderr, "can not read pre-computed parameters file\n");
	}
	else {
		auto tmp = std::string(filename);
		auto p1 = tmp.find_last_of('/') + 1;
		auto p2 = tmp.find_last_of('.');
		auto mat_name = tmp.substr(p1, p2 - p1);
		std::string params;
		std::regex re(mat_name +
					  R"(\s([0-9]*\.?[0-9]+)\s\((\S+),(\S+),(\S+)\)\s\((\S+),(\S+),(\S+)\)\s\((\S+),(\S+),(\S+)\)\s(\d+)(\**))");
		while (std::getline(paramsfs, params)) {
			std::smatch matches;
			if (std::regex_match(params, matches, re)) {
				if (cmbParamsPack == nullptr) {
					cmbParamsPack = new MERLParametersPack();
				}
				cmbParamsPack->ggx = std::atof(matches[1].str().c_str());
				cmbParamsPack->diffuseAlbedo[0] = std::atof(matches[2].str().c_str());
				cmbParamsPack->diffuseAlbedo[1] = std::atof(matches[3].str().c_str());
				cmbParamsPack->diffuseAlbedo[2] = std::atof(matches[4].str().c_str());
				cmbParamsPack->specularAlbedo[0] = std::atof(matches[5].str().c_str());
				cmbParamsPack->specularAlbedo[1] = std::atof(matches[6].str().c_str());
				cmbParamsPack->specularAlbedo[2] = std::atof(matches[7].str().c_str());
				cmbParamsPack->optimalThreshold[0] = std::atof(matches[8].str().c_str());
				cmbParamsPack->optimalThreshold[1] = std::atof(matches[9].str().c_str());
				cmbParamsPack->optimalThreshold[2] = std::atof(matches[10].str().c_str());
				cmbParamsPack->cluster = std::atoi(matches[11].str().c_str());

				// reset specular parameters
				getFloatParameter(3)->defaultVal = cmbParamsPack->specularAlbedo[0];
				getFloatParameter(4)->defaultVal = cmbParamsPack->specularAlbedo[1];
				getFloatParameter(5)->defaultVal = cmbParamsPack->specularAlbedo[2];
			
				//getFloatParameter(3)->currentVal = getFloatParameter(3)->defaultVal;
				//getFloatParameter(4)->currentVal = getFloatParameter(4)->defaultVal;
				//getFloatParameter(5)->currentVal = getFloatParameter(5)->defaultVal;
			}
		}
	}

	cmbEnabled = true;
	return true;
}


void BRDFMeasuredMERL::enableCmb(bool enable) {
	if (cmbParamsPack == nullptr) return;
	cmbEnabled = enable;
	if (enable) {
		getFloatParameter(3)->defaultVal = cmbParamsPack->specularAlbedo[0];
		getFloatParameter(4)->defaultVal = cmbParamsPack->specularAlbedo[1];
		getFloatParameter(5)->defaultVal = cmbParamsPack->specularAlbedo[2];
	}
	else {
		// restore
		getFloatParameter(3)->defaultVal = paramsPack->specularAlbedo[0];
		getFloatParameter(4)->defaultVal = paramsPack->specularAlbedo[1];
		getFloatParameter(5)->defaultVal = paramsPack->specularAlbedo[2];
		getFloatParameter(3)->currentVal = getFloatParameter(3)->defaultVal;
		getFloatParameter(4)->currentVal = getFloatParameter(4)->defaultVal;
		getFloatParameter(5)->currentVal = getFloatParameter(5)->defaultVal;
	}
}


void BRDFMeasuredMERL::initGL()
{
    if( initializedGL )
        return;
    
    
    // create buffer object
    glGenBuffers(1, &tbo);
    glBindBuffer(GL_TEXTURE_BUFFER_EXT, tbo);
    
    
    
    // initialize buffer object
    unsigned int numBytes = numBRDFSamples * 3 * sizeof(float);
    //printf( "size = %d bytes (%f megs)\n", numBytes, float(numBytes) / 1048576.0f );
    glBufferData( GL_TEXTURE_BUFFER_EXT, numBytes, 0, GL_STATIC_DRAW );
    
    //tex
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_BUFFER_EXT, tex);
    glTexBufferEXT(GL_TEXTURE_BUFFER_EXT, GL_INTENSITY32F_ARB, tbo);
    glBindBuffer(GL_TEXTURE_BUFFER_EXT, 0);
    
    
    
    glBindBuffer(GL_TEXTURE_BUFFER_EXT, tbo);
    float* p = (float*)glMapBuffer( GL_TEXTURE_BUFFER_EXT, GL_WRITE_ONLY );
    
    
    memcpy( p, brdfData, numBytes );
    glUnmapBuffer(GL_TEXTURE_BUFFER_EXT);
    glBindBuffer(GL_TEXTURE_BUFFER_EXT, 0);
    
    // keep origin data for later use
    //delete[] brdfData;
    //brdfData = NULL;
    
    initializedGL = true;
}


void BRDFMeasuredMERL::adjustShaderPreRender( DGLShader* shader )
{
    shader->setUniformTexture( "measuredData", tex, GL_TEXTURE_BUFFER_EXT );

	shader->setUniformInt("enableCmb", cmbEnabled);
	
	if (cmbTEX) {
		shader->setUniformTexture("cmbData", cmbTEX, GL_TEXTURE_BUFFER_EXT);
	}

	if (cmbParamsPack != nullptr) {
		shader->setUniformFloat("cmbThreshold", cmbParamsPack->optimalThreshold[0],
								cmbParamsPack->optimalThreshold[1], cmbParamsPack->optimalThreshold[2]);
	}

	if (paramsPack != nullptr) {
		shader->setUniformFloat("optimalThreshold", paramsPack->optimalThreshold[0],
								paramsPack->optimalThreshold[1], paramsPack->optimalThreshold[2]);

		shader->setUniformFloat("diffuseAlbedoRatio",
								getFloatParameter(0)->currentVal / getFloatParameter(0)->defaultVal,
								getFloatParameter(1)->currentVal / getFloatParameter(1)->defaultVal,
								getFloatParameter(2)->currentVal / getFloatParameter(2)->defaultVal);

		if (getFloatParameter(3)->defaultVal == 0) {
			shader->setUniformFloat("specularAlbedoRatio", 0, 0, 0);
		}
		else {
			shader->setUniformFloat("specularAlbedoRatio",
									getFloatParameter(3)->currentVal / getFloatParameter(3)->defaultVal,
									getFloatParameter(4)->currentVal / getFloatParameter(4)->defaultVal,
									getFloatParameter(5)->currentVal / getFloatParameter(5)->defaultVal);
		}

	}

    BRDFBase::adjustShaderPreRender( shader );
}


