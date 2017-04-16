#ifndef MERL_PARAMETERS_PACK_H
#define MERL_PARAMETERS_PACK_H

#include <string>

struct MERLParametersPack {

	float ggx;
	float diffuseAlbedo[3];
	float specularAlbedo[3];
	float optimalThreshold[3];
	int cluster;

};

#endif
