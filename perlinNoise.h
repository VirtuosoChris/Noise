#ifndef PERLIN_NOISE_H_INCLUDED
#define PERLIN_NOISE_H_INCLUDED

/*
MIT License

Copyright(c) 2012 Christopher Pugh

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and /or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

///\todo replace srand(time(0)) and other random calls with better uniform random generation

#include <Eigen/Geometry>
#include <ctime>

#define PERLIN_TABLE_SIZE 0x100 //256 entries in table

template <class T>
inline T s_curve(T t)
{
    const float tSquared = t*t;
    return 3.0 * tSquared - 2.0 * t * tSquared;
}

template<class T>
inline T lerp(T a, T b, T t)
{
    return (1.0 - t) * a + t*b;
}

template<class T>
inline T s_interpolate(T a, T b, T t)
{
    return lerp(a, b, s_curve(t));
}

inline float randomFloatInRange(const float& min, const float& max)
{
    return (((float)rand() / (float)RAND_MAX) * (max-min)) + min;
}

template<int DIMENSIONS>
class PerlinGenerator
{
public:

    std::vector<Eigen::Vector2d> gradients;
    std::vector<unsigned int> permutations;

    PerlinGenerator()
    {
        gradients.reserve(PERLIN_TABLE_SIZE);
        permutations.reserve(PERLIN_TABLE_SIZE);

        //create our permutations table and stuff
        for(unsigned int i = 0; i < PERLIN_TABLE_SIZE; i++)
        {
            permutations.push_back( rand() % PERLIN_TABLE_SIZE);

            gradients[i] = Eigen::Vector2d(randomFloatInRange(-1.0f, 1.0f), randomFloatInRange(-1.0f, 1.0f));
            gradients[i].normalize();
        }
    }

    double calculateNoise(Eigen::Vector2d atPos, double wavelength)
    {
        const double frequency = 1.0 / wavelength;
        const double invRoot2 = 1.0 / sqrt(2.0); //normalizing factor since vector in unit square has max length root 2

        Eigen::Vector2d dArrayIndices = atPos * frequency;

        Eigen::Vector2d baseIndices = Eigen::Vector2d( std::floor(dArrayIndices[0]), std::floor(dArrayIndices[1]));

        Eigen::Vector2d lerpParams = dArrayIndices - baseIndices;

        Eigen::Vector2i tableIndices =  baseIndices.cast<int>(); //cast to int

        tableIndices[0] &= 0xff;
        tableIndices[1] &= 0xff;

        Eigen::Vector2i incrementedIndices = tableIndices;
        incrementedIndices[0]++; incrementedIndices[1]++;

        //mod with 256.  The and automagically takes care of getting rid of the negative table indices we might otherwise get.  Thanks Ken Perlin
        //and two's complement!
        incrementedIndices[0] &= 0xff;
        incrementedIndices[1] &= 0xff;

        const int topLeftIndex = (tableIndices[0] + permutations[incrementedIndices[1]]) & 0xff;
        const int topRightIndex = (incrementedIndices[0] + permutations[incrementedIndices[1]]) & 0xff;
        const int bottomLeftIndex = (tableIndices[0] + permutations[tableIndices[1]])&0xff;
        const int bottomRightIndex = (incrementedIndices[0] + permutations[tableIndices[1]])&0xff;

        const Eigen::Vector2d& upperLeftGrad = gradients[topLeftIndex];
        const Eigen::Vector2d& upperRightGrad = gradients[topRightIndex];
        const Eigen::Vector2d& lowerLeftGrad = gradients[bottomLeftIndex];
        const Eigen::Vector2d& lowerRightGrad = gradients[bottomRightIndex];

        ///\todo : scale, setting shit right
        Eigen::Vector2d tempVec4 = lerpParams;//lowerLeftCorner-atPos;
        Eigen::Vector2d tempVec1 = Eigen::Vector2d(lerpParams[0], -(1.0 - lerpParams[1]));//upperLeftCorner-atPos;
        Eigen::Vector2d tempVec2 = Eigen::Vector2d(- (1.0 - lerpParams[0]), -(1.0 - lerpParams[1]));//upperRightCorner-atPos;
        Eigen::Vector2d tempVec3 = Eigen::Vector2d(- (1.0 - lerpParams[0]), lerpParams[1]);  //lowerRightCorner-atPos;

        double ulVal = (tempVec1).dot(upperLeftGrad);
        double urVal = (tempVec2).dot(upperRightGrad);
        double lrVal = (tempVec3).dot(lowerRightGrad);
        double llVal = (tempVec4).dot(lowerLeftGrad);

        double lerpXTop =  s_interpolate(ulVal, urVal, lerpParams[0]);
        double lerpXBottom = s_interpolate(llVal, lrVal, lerpParams[0]);

        double result = s_interpolate(lerpXBottom, lerpXTop, lerpParams[1]);

        ///\todo this should never have to be done!
        result = std::min(1.0, result);
        result = std::max(-1.0,result);

        return result;
    }
};

inline void generatePerlinNoiseSlice(unsigned char* image, unsigned int width, unsigned int height, float amplitude, float wavelength)
{
    //get gradient vectors
    const float frequency = 1.0 / wavelength;
    const unsigned int gradSize = 256;

    Eigen::Vector2f gradientArrays[gradSize][gradSize];
    srand(time(0));

    //generate gradient vectors
    for(int i = 0; i <gradSize; i++)
    {
        for(int j = 0; j < gradSize; j++)
        {
            gradientArrays[i][j] = Eigen::Vector2f(randomFloatInRange(-1.0f, 1.0f), randomFloatInRange(-1.0f, 1.0f));
            gradientArrays[i][j].normalize();
        }
    }

    for(unsigned int h = 0; h < height; h++)
    {
        for(unsigned int w = 0; w < width; w++)
        {
            float xOff = float(w);
            float yOff = float(h);

            Eigen::Vector2f atPos = Eigen::Vector2f(xOff, yOff);

            float tmpA = xOff * frequency;
            float tmpB = yOff * frequency;

            int leftIndex =  (int)std::floor(tmpA);
            int rightIndex = leftIndex+1;
            int bottomIndex = (int)std::floor(tmpB);
            int topIndex = bottomIndex+1;

            ///todo bottom wraparound
            int leftIndexMod = leftIndex % gradSize;
            int rightIndexMod = rightIndex% gradSize;
            int bottomIndexMod = bottomIndex% gradSize;
            int topIndexMod = topIndex% gradSize;

            Eigen::Vector2f& upperLeftGrad = gradientArrays[topIndexMod][leftIndexMod];
            Eigen::Vector2f& upperRightGrad = gradientArrays[topIndexMod][rightIndexMod];
            Eigen::Vector2f& lowerLeftGrad = gradientArrays[bottomIndexMod][leftIndexMod];
            Eigen::Vector2f& lowerRightGrad = gradientArrays[bottomIndexMod][rightIndexMod];

            float leftPos = leftIndex * wavelength;
            float bottomPos = bottomIndex * wavelength;

            Eigen::Vector2f upperLeftCorner = Eigen::Vector2f(leftPos , bottomPos + wavelength);
            Eigen::Vector2f lowerLeftCorner = Eigen::Vector2f(leftPos, bottomPos);
            Eigen::Vector2f upperRightCorner = Eigen::Vector2f(leftPos + wavelength, bottomPos + wavelength);
            Eigen::Vector2f lowerRightCorner = Eigen::Vector2f(leftPos + wavelength, bottomPos);

            Eigen::Vector2f tempVec1 = upperLeftCorner-atPos;
            Eigen::Vector2f tempVec2 = upperRightCorner-atPos;
            Eigen::Vector2f tempVec3 = lowerRightCorner-atPos;
            Eigen::Vector2f tempVec4 = lowerLeftCorner-atPos;

            const float invRoot2 = 1.0f / sqrt(2.0);

            tempVec1 *= frequency * invRoot2;
            tempVec2 *= frequency * invRoot2;
            tempVec3 *= frequency * invRoot2;
            tempVec4 *= frequency * invRoot2;

            float ulVal = (tempVec1).dot(upperLeftGrad);
            float urVal = (tempVec2).dot(upperRightGrad);
            float lrVal = (tempVec3).dot(lowerRightGrad);
            float llVal = (tempVec4).dot(lowerLeftGrad);

            float xLerpVal = (atPos[0] - upperLeftCorner[0]) / (upperRightCorner[0] - upperLeftCorner[0]);

            float lerpXTop =  s_interpolate(ulVal, urVal, xLerpVal);
            float lerpXBottom = s_interpolate(llVal, lrVal, xLerpVal);

            float yLerpVal = (atPos[1] - lowerLeftCorner[1]) / (upperLeftCorner[1] - lowerLeftCorner[1]);

            float result = s_interpolate(lerpXBottom, lerpXTop, yLerpVal);

            result = std::min(1.0f, result);
            result = std::max(-1.0f,result);

            result = ((result + 1.0) * .5);

            image[h * width + w] = (unsigned char)(result * 255.0f);
        }
    }
}

#endif
