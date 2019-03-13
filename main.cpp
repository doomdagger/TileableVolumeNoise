

#include <iostream>
#include <fstream>
#include <cinttypes>
#include <time.h>
#include <math.h>
#include <string>

#include "./TileableVolumeNoise.h"
#include "./libtarga.h"

#include <ppl.h>
using namespace concurrency;

void writeTGA(const char* fileName, int width, int height, /*const*/ unsigned char* data)
{
	if (!tga_write_raw(fileName, width, height, data, TGA_TRUECOLOR_32))
	{
		printf("Failed to write image!\n");
		printf(tga_error_string(tga_get_last_error()));
	}
}

void writeKTX(const std::string &fileName, int origWidth, int origHeight, int origDepth, unsigned char *data, bool tile3D=false)
{
	int texDepth = origDepth;
	int texWidth = origWidth;
	int texHeight = origHeight;
	int channels = 4;
	bool tileEnabled = tile3D && texDepth > 1;
	if (tileEnabled)
	{
		texWidth *= texDepth;
		texDepth = 1;
	}
	int texSize = texWidth * texHeight * texDepth;

	byte header[] = {
		(char)0xAB, (char)0x4B, (char)0x54, (char)0x58, // first four bytes of Byte[12] identifier
		(char)0x20, (char)0x31, (char)0x31, (char)0xBB, // next four bytes of Byte[12] identifier
		(char)0x0D, (char)0x0A, (char)0x1A, (char)0x0A, // final four bytes of Byte[12] identifier
		(char)0x01, (char)0x02, (char)0x03, (char)0x04, // Byte[4] endianess (Big endian in this case)
	};

	std::ofstream fs(fileName, std::ios::binary);
	fs.write(header, 16);

	uint32_t glType = 0x1401; // unsigned byte
	uint32_t glTypeSize = 1; // 1 byte
	uint32_t glFormat = 0x1908; // RGBA
	uint32_t glInterformat = 0x8D7C; // RGBA8UI
	uint32_t glBaseInternalFormat = 0x1908; // RGBA
	uint32_t width = (uint32_t)texWidth;
	uint32_t height = (uint32_t)texHeight;
	uint32_t depth = (uint32_t)(texDepth == 1 ? 0 : texDepth);
	uint32_t numOfArrElem = 0;
	uint32_t numOfFace = 1;
	uint32_t numOfMip = 1;
	uint32_t bytesOfKeyVal = 0;

	fs.write(reinterpret_cast<const char*>(&glType), sizeof uint32_t);
	fs.write(reinterpret_cast<const char*>(&glTypeSize), sizeof uint32_t);
	fs.write(reinterpret_cast<const char*>(&glFormat), sizeof uint32_t);
	fs.write(reinterpret_cast<const char*>(&glInterformat), sizeof uint32_t);
	fs.write(reinterpret_cast<const char*>(&glBaseInternalFormat), sizeof uint32_t);
	fs.write(reinterpret_cast<const char*>(&width), sizeof uint32_t);
	fs.write(reinterpret_cast<const char*>(&height), sizeof uint32_t);
	fs.write(reinterpret_cast<const char*>(&depth), sizeof uint32_t);
	fs.write(reinterpret_cast<const char*>(&numOfArrElem), sizeof uint32_t);
	fs.write(reinterpret_cast<const char*>(&numOfFace), sizeof uint32_t);
	fs.write(reinterpret_cast<const char*>(&numOfMip), sizeof uint32_t);
	fs.write(reinterpret_cast<const char*>(&bytesOfKeyVal), sizeof uint32_t);

	uint32_t imageSize = (uint32_t)(texSize * channels * glTypeSize);
	fs.write(reinterpret_cast<const char*>(&imageSize), sizeof uint32_t);

	if (tileEnabled)
	{
		for (int j = 0; j < origHeight; j++)
			for (int k = 0; k < origDepth; k++)
				for (int i = 0; i < origWidth; i++)
				{
					int startIndex = k * origWidth * origHeight * channels + j * origWidth * channels + i * channels;
					fs.write(reinterpret_cast<const char*>(&data[startIndex]), 4);
				}
	}
	else
	{
		int dataSize = texSize * channels;
		fs.write(reinterpret_cast<const char*>(data), texSize * channels);
	}

	fs.close();
}

// the remap function used in the shaders as described in Gpu Pro 7. It must match when using pre packed textures
float remap(float originalValue, float originalMin, float originalMax, float newMin, float newMax)
{
	return newMin + (((originalValue - originalMin) / (originalMax - originalMin)) * (newMax - newMin));
}

int main (int argc, char *argv[])
{   
	//
	// Exemple of tileable Perlin noise texture generation
	//

	/*
	unsigned int gPerlinNoiseTextureSize = 32;
	unsigned char* perlinNoiseTexels = (unsigned char*)malloc(gPerlinNoiseTextureSize*gPerlinNoiseTextureSize*gPerlinNoiseTextureSize * sizeof(unsigned char));

	// Generate Perlin noise source
	const glm::vec3 normFactPerlin = glm::vec3(1.0f / float(gPerlinNoiseTextureSize));
	parallel_for(int(0), int(gPerlinNoiseTextureSize), [&](int s) //for (int s=0; s<giPerlinNoiseTextureSize; s++)
	{
		for (int t = 0; t<gPerlinNoiseTextureSize; t++)
		{
			for (int r = 0; r<gPerlinNoiseTextureSize; r++)
			{
				glm::vec3 coord = glm::vec3(s, t, r) * normFactPerlin;

				const int octaveCount = 1;
				const float frequency = 8;
				float noise = Tileable3dNoise::PerlinNoise(coord, frequency, octaveCount);
				noise *= 255.0f;

				int addr = r*gPerlinNoiseTextureSize*gPerlinNoiseTextureSize + t*gPerlinNoiseTextureSize + s;
				perlinNoiseTexels[addr] = unsigned char(noise);

			}
		}
	}
	); // end parallel_for
	free(perlinNoiseTexels);

	//
	// Exemple of tileable Worley noise texture generation
	//
	unsigned int gWorleyNoiseTextureSize = 32;
	unsigned char* worleyNoiseTexels = (unsigned char*)malloc(gWorleyNoiseTextureSize*gWorleyNoiseTextureSize*gWorleyNoiseTextureSize * sizeof(unsigned char));

	const glm::vec3 normFactWorley = glm::vec3(1.0f / float(gWorleyNoiseTextureSize));
	parallel_for(int(0), int(gWorleyNoiseTextureSize), [&](int s) //for (int s = 0; s<giWorleyNoiseTextureSize; s++)
	{
		for (int t = 0; t<gWorleyNoiseTextureSize; t++)
		{
			for (int r = 0; r<gWorleyNoiseTextureSize; r++)
			{
				glm::vec3 coord = glm::vec3(s, t, r) * normFactPerlin;

				const float cellCount = 3;
				float noise = 1.0 - Tileable3dNoise::WorleyNoise(coord, cellCount);
				noise *= 255.0f;

				int addr = r*gWorleyNoiseTextureSize*gWorleyNoiseTextureSize + t*gWorleyNoiseTextureSize + s;
				worleyNoiseTexels[addr] = unsigned char(noise);
			}
		}
	}
	); // end parallel_for
	free(worleyNoiseTexels);
	*/



	//
	// Generate cloud shape and erosion texture similarly GPU Pro 7 chapter II-4
	//

	// Frequence multiplicator. No boudary check etc. but fine for this small tool.
	const float frequenceMul[6] = { 2.0f,8.0f,14.0f,20.0f,26.0f,32.0f };	// special weight for perling worley

														// Cloud base shape (will be used to generate PerlingWorley noise in he shader)
														// Note: all channels could be combined once here to reduce memory bandwith requirements.
	int cloudBaseShapeTextureSize = 128;				// !!! If this is reduce, you hsould also reduce the number of frequency in the fmb noise  !!!
	int cloudBaseShapeRowBytes = cloudBaseShapeTextureSize * sizeof(unsigned char) * 4;
	int cloudBaseShapeSliceBytes = cloudBaseShapeRowBytes * cloudBaseShapeTextureSize;
	int cloudBaseShapeVolumeBytes = cloudBaseShapeSliceBytes * cloudBaseShapeTextureSize;
	unsigned char* cloudBaseShapeTexels = (unsigned char*)malloc(cloudBaseShapeVolumeBytes);
	//unsigned char* cloudBaseShapeTexelsPacked = (unsigned char*)malloc(cloudBaseShapeVolumeBytes);
	parallel_for(int(0), int(cloudBaseShapeTextureSize), [&](int s) //for (int s = 0; s<gCloudBaseShapeTextureSize; s++)
	{
		const glm::vec3 normFact = glm::vec3(1.0f / float(cloudBaseShapeTextureSize));
		for (int t = 0; t<cloudBaseShapeTextureSize; t++)
		{
			for (int r = 0; r<cloudBaseShapeTextureSize; r++)
			{
				glm::vec3 coord = glm::vec3(s, t, r) * normFact;

				// Perlin FBM noise
				const int octaveCount = 3;
				const float frequency = 8.0f;
				float perlinNoise = Tileable3dNoise::PerlinNoise(coord, frequency, octaveCount);

				float PerlinWorleyNoise = 0.0f;
				{
					const float cellCount = 4;
					const float worleyNoise0 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * frequenceMul[0]));
					const float worleyNoise1 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * frequenceMul[1]));
					const float worleyNoise2 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * frequenceMul[2]));
					const float worleyNoise3 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * frequenceMul[3]));
					const float worleyNoise4 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * frequenceMul[4]));
					const float worleyNoise5 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * frequenceMul[5]));	// half the frequency of texel, we should not go further (with cellCount = 32 and texture size = 64)

					float worleyFBM = worleyNoise0*0.625f + worleyNoise1*0.25f + worleyNoise2*0.125f;

					// Perlin Worley is based on description in GPU Pro 7: Real Time Volumetric Cloudscapes.
					// However it is not clear the text and the image are matching: images does not seem to match what the result  from the description in text would give.
					// Also there are a lot of fudge factor in the code, e.g. *0.2, so it is really up to you to fine the formula you like.
					PerlinWorleyNoise = remap(worleyFBM, 0.0, 1.0, 0.0, perlinNoise);	// Matches better what figure 4.7 (not the following up text description p.101). Maps worley between newMin as 0 and 
					//PerlinWorleyNoise = remap(perlinNoise, 0.0f, 1.0f, worleyFBM, 1.0f);	// mapping perlin noise in between worley as minimum and 1.0 as maximum (as described in text of p.101 of GPU Pro 7) 
				}

				const float cellCount = 4;
				float worleyNoise0 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * 1));
				float worleyNoise1 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * 2));
				float worleyNoise2 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * 4));
				float worleyNoise3 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * 8));
				float worleyNoise4 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * 16));
				//float worleyNoise5 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * 32));	//cellCount=2 -> half the frequency of texel, we should not go further (with cellCount = 32 and texture size = 64)

				// Three frequency of Worley FBM noise
				float worleyFBM0 = worleyNoise1*0.625f + worleyNoise2*0.25f + worleyNoise3*0.125f;
				float worleyFBM1 = worleyNoise2*0.625f + worleyNoise3*0.25f + worleyNoise4*0.125f;
				//float worleyFBM2 = worleyNoise3*0.625f + worleyNoise4*0.25f + worleyNoise5*0.125f;
				float worleyFBM2 = worleyNoise3*0.75f + worleyNoise4*0.25f; // cellCount=4 -> worleyNoise5 is just noise due to sampling frequency=texel frequency. So only take into account 2 frequencies for FBM

				int addr = r*cloudBaseShapeTextureSize*cloudBaseShapeTextureSize + t*cloudBaseShapeTextureSize + s;

				addr *= 4;
				cloudBaseShapeTexels[addr] = unsigned char(255.0f*PerlinWorleyNoise);
				cloudBaseShapeTexels[addr + 1] = unsigned char(255.0f*worleyFBM0);
				cloudBaseShapeTexels[addr + 2] = unsigned char(255.0f*worleyFBM1);
				cloudBaseShapeTexels[addr + 3] = unsigned char(255.0f*worleyFBM2);

				//float value = 0.0;
				//{
				//	// pack the channels for direct usage in shader
				//	float lowFreqFBM = worleyFBM0*0.625f + worleyFBM1*0.25f + worleyFBM2*0.125f;
				//	float baseCloud = PerlinWorleyNoise;
				//	value = remap(baseCloud, -(1.0f - lowFreqFBM), 1.0f, 0.0f, 1.0f);
				//	// Saturate
				//	value = std::fminf(value, 1.0f);
				//	value = std::fmaxf(value, 0.0f);
				//}
				//cloudBaseShapeTexelsPacked[addr] = unsigned char(255.0f*value);
				//cloudBaseShapeTexelsPacked[addr + 1] = unsigned char(255.0f*value);
				//cloudBaseShapeTexelsPacked[addr + 2] = unsigned char(255.0f*value);
				//cloudBaseShapeTexelsPacked[addr + 3] = unsigned char(255.0f);
			}
		}
	}
	); // end parallel_for
	{
		int width = cloudBaseShapeTextureSize*cloudBaseShapeTextureSize;
		int height = cloudBaseShapeTextureSize;
		writeKTX("noiseShape.ktx", cloudBaseShapeTextureSize, cloudBaseShapeTextureSize, cloudBaseShapeTextureSize, cloudBaseShapeTexels);
		//writeTGA("noiseShape.tga",       width, height, cloudBaseShapeTexels);
		//writeTGA("noiseShapePacked.tga", width, height, cloudBaseShapeTexelsPacked);
	}






	// Detail texture behing different frequency of Worley noise
	// Note: all channels could be combined once here to reduce memory bandwith requirements.
	int cloudErosionTextureSize = 32;
	int cloudErosionRowBytes = cloudErosionTextureSize * sizeof(unsigned char) * 4;
	int cloudErosionSliceBytes = cloudErosionRowBytes * cloudErosionTextureSize;
	int cloudErosionVolumeBytes = cloudErosionSliceBytes * cloudErosionTextureSize;
	unsigned char* cloudErosionTexels = (unsigned char*)malloc(cloudErosionVolumeBytes);
	//unsigned char* cloudErosionTexelsPacked = (unsigned char*)malloc(cloudErosionVolumeBytes);
	parallel_for(int(0), int(cloudErosionTextureSize), [&](int s) //for (int s = 0; s<gCloudErosionTextureSize; s++)
	{
		const glm::vec3 normFact = glm::vec3(1.0f / float(cloudErosionTextureSize));
		for (int t = 0; t<cloudErosionTextureSize; t++)
		{
			for (int r = 0; r<cloudErosionTextureSize; r++)
			{
				glm::vec3 coord = glm::vec3(s, t, r) * normFact;

#if 1
				// 3 octaves
				const float cellCount = 2;
				float worleyNoise0 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * 1));
				float worleyNoise1 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * 2));
				float worleyNoise2 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * 4));
				float worleyNoise3 = (1.0f - Tileable3dNoise::WorleyNoise(coord, cellCount * 8));
				float worleyFBM0 = worleyNoise0*0.625f + worleyNoise1*0.25f + worleyNoise2*0.125f;
				float worleyFBM1 = worleyNoise1*0.625f + worleyNoise2*0.25f + worleyNoise3*0.125f;
				float worleyFBM2 = worleyNoise2*0.75f + worleyNoise3*0.25f; // cellCount=4 -> worleyNoise4 is just noise due to sampling frequency=texel freque. So only take into account 2 frequencies for FBM
#else
				// 2 octaves
				float worleyNoise0 = (1.0f - Tileable3dNoise::WorleyNoise(coord, 4));
				float worleyNoise1 = (1.0f - Tileable3dNoise::WorleyNoise(coord, 7));
				float worleyNoise2 = (1.0f - Tileable3dNoise::WorleyNoise(coord, 10));
				float worleyNoise3 = (1.0f - Tileable3dNoise::WorleyNoise(coord, 13));
				float worleyFBM0 = worleyNoise0*0.75f + worleyNoise1*0.25f;
				float worleyFBM1 = worleyNoise1*0.75f + worleyNoise2*0.25f;
				float worleyFBM2 = worleyNoise2*0.75f + worleyNoise3*0.25f;
#endif

				int addr = r*cloudErosionTextureSize*cloudErosionTextureSize + t*cloudErosionTextureSize + s;
				addr *= 4;
				cloudErosionTexels[addr] = unsigned char(255.0f*worleyFBM0);
				cloudErosionTexels[addr + 1] = unsigned char(255.0f*worleyFBM1);
				cloudErosionTexels[addr + 2] = unsigned char(255.0f*worleyFBM2);
				cloudErosionTexels[addr + 3] = unsigned char(255.0f);

				//float value = 0.0;
				//{
				//	value = worleyFBM0*0.625f + worleyFBM1*0.25f + worleyFBM2*0.125f;
				//}
				//cloudErosionTexelsPacked[addr] = unsigned char(255.0f * value);
				//cloudErosionTexelsPacked[addr + 1] = unsigned char(255.0f * value);
				//cloudErosionTexelsPacked[addr + 2] = unsigned char(255.0f * value);
				//cloudErosionTexelsPacked[addr + 3] = unsigned char(255.0f);
			}
		}
	}
	); // end parallel_for
	{
		int width = cloudErosionTextureSize*cloudErosionTextureSize;
		int height = cloudErosionTextureSize;
		writeKTX("noiseErosion.ktx", cloudErosionTextureSize, cloudErosionTextureSize, cloudErosionTextureSize, cloudErosionTexels);
		//writeTGA("noiseErosion.tga",       width, height, cloudErosionTexels);
		//writeTGA("noiseErosionPacked.tga", width, height, cloudErosionTexelsPacked);
	}

#if 0
	// Debug tileability using a 3x3 tile of the same slice, see if edges appears.

	auto debugPrintTileability = [&](auto addrSrc2, auto addrDst2, const char* debugStr) 
	{
		for (int r = 0; r < cloudBaseShapeTextureSize; r += cloudBaseShapeTextureSize / 8)
		{
			unsigned char* debugImg = (unsigned char*)malloc(9 * cloudBaseShapeSliceBytes);
			for (int i = 0; i < cloudBaseShapeTextureSize; i++)
			{
				int t = i;

				// copy the row into the 9 debug texture tiles
				auto addrSrc = [&]()
				{
					return addrSrc2(r, t, cloudBaseShapeSliceBytes, cloudBaseShapeRowBytes);
				};
				auto addrDst = [&](auto x, auto y)
				{
					return addrDst2(x, y, t, cloudBaseShapeSliceBytes, cloudBaseShapeRowBytes);
				};

				memcpy(debugImg + addrDst(0, 0), cloudBaseShapeTexelsPacked + addrSrc(), cloudBaseShapeRowBytes);
				memcpy(debugImg + addrDst(1, 0), cloudBaseShapeTexelsPacked + addrSrc(), cloudBaseShapeRowBytes);
				memcpy(debugImg + addrDst(2, 0), cloudBaseShapeTexelsPacked + addrSrc(), cloudBaseShapeRowBytes);
				memcpy(debugImg + addrDst(0, 1), cloudBaseShapeTexelsPacked + addrSrc(), cloudBaseShapeRowBytes);
				memcpy(debugImg + addrDst(1, 1), cloudBaseShapeTexelsPacked + addrSrc(), cloudBaseShapeRowBytes);
				memcpy(debugImg + addrDst(2, 1), cloudBaseShapeTexelsPacked + addrSrc(), cloudBaseShapeRowBytes);
				memcpy(debugImg + addrDst(0, 2), cloudBaseShapeTexelsPacked + addrSrc(), cloudBaseShapeRowBytes);
				memcpy(debugImg + addrDst(1, 2), cloudBaseShapeTexelsPacked + addrSrc(), cloudBaseShapeRowBytes);
				memcpy(debugImg + addrDst(2, 2), cloudBaseShapeTexelsPacked + addrSrc(), cloudBaseShapeRowBytes);


			}
			char fileName[256];
			sprintf_s(fileName, 256, "debugBase%s%i.tga", debugStr, r);
			writeTGA(fileName, cloudBaseShapeTextureSize * 3, cloudBaseShapeTextureSize * 3, debugImg);
		}

		for (int r = 0; r < cloudErosionTextureSize; r += cloudErosionTextureSize / 8)
		{
			unsigned char* debugImg = (unsigned char*)malloc(9*cloudErosionSliceBytes);
			for (int i = 0; i < cloudErosionTextureSize; i++)
			{
				int t = i;

				auto addrSrc = [&]()
				{
					return addrSrc2(r, t, cloudErosionSliceBytes, cloudErosionRowBytes);
				};
				auto addrDst = [&](auto x, auto y)
				{
					return addrDst2(x, y, t, cloudErosionSliceBytes, cloudErosionRowBytes);
				};

				memcpy(debugImg + addrDst(0,0), cloudErosionTexelsPacked + addrSrc(), cloudErosionRowBytes);
				memcpy(debugImg + addrDst(1,0), cloudErosionTexelsPacked + addrSrc(), cloudErosionRowBytes);
				memcpy(debugImg + addrDst(2,0), cloudErosionTexelsPacked + addrSrc(), cloudErosionRowBytes);
				memcpy(debugImg + addrDst(0,1), cloudErosionTexelsPacked + addrSrc(), cloudErosionRowBytes);
				memcpy(debugImg + addrDst(1,1), cloudErosionTexelsPacked + addrSrc(), cloudErosionRowBytes);
				memcpy(debugImg + addrDst(2,1), cloudErosionTexelsPacked + addrSrc(), cloudErosionRowBytes);
				memcpy(debugImg + addrDst(0,2), cloudErosionTexelsPacked + addrSrc(), cloudErosionRowBytes);
				memcpy(debugImg + addrDst(1,2), cloudErosionTexelsPacked + addrSrc(), cloudErosionRowBytes);
				memcpy(debugImg + addrDst(2,2), cloudErosionTexelsPacked + addrSrc(), cloudErosionRowBytes);
			}
			char fileName[256];
			sprintf_s(fileName, 256, "debugErosion%s%i.tga", debugStr, r);
			writeTGA(fileName, cloudErosionTextureSize*3, cloudErosionTextureSize*3, debugImg);
		}
	};

	auto addrDst = [](auto x, auto y, auto t, auto sliceBytes, auto rowBytes)
	{
		return x * rowBytes + y * sliceBytes * 3 + t * rowBytes * 3;
	};
	auto addrSrcXY = [](auto r, auto t, auto sliceBytes, auto rowBytes)
	{
		return r*sliceBytes + t*rowBytes;
	};
	auto addrSrcXZ = [](auto r, auto t, auto sliceBytes, auto rowBytes)
	{
		return t*sliceBytes + r*rowBytes;
	};

	debugPrintTileability(addrSrcXY, addrDst, "XY");
	debugPrintTileability(addrSrcXZ, addrDst, "XZ");
#endif

	free(cloudBaseShapeTexels);
	//free(cloudErosionTexelsPacked);
	free(cloudErosionTexels);
	//free(cloudErosionTexelsPacked);

    return 0;
}

