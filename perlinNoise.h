#pragma once

/*
MIT License

Copyright(c) 2012-2024 Christopher Pugh

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

/*

todos:
havent tested this heavily refactored version yet
arbitrary dimensionality
shaders (fixed 2d)
shaders (generated?)
other noise types
threading
port sample proc gen examples
*/

#include <Eigen/Geometry>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <concepts>

template <class T>
inline T s_curve(T t)
{
    const T tSquared = t*t;
    return T(3.0) * tSquared - T(2.0) * t * tSquared;
}

template<class T>
inline T lerp(T a, T b, T t)
{
    return (T(1.0) - t) * a + t*b;
}

template<class T>
inline T s_interpolate(T a, T b, T t)
{
    return lerp(a, b, s_curve(t));
}

template<typename Functor, typename RealType>
concept ValidFunctor = requires(Functor f, RealType value)
{
    { f(value) } -> std::same_as<RealType>;
};

template <typename RealType = double>
requires std::is_floating_point_v<RealType>
class PerlinGenerator
{
private:

    static constexpr std::size_t PERLIN_TABLE_SIZE = 0x100;

    using Vector2 = Eigen::Matrix<RealType, 2, 1>;

    std::vector<Vector2> gradients;
    std::vector<unsigned int> permutations;

    void initialize(std::default_random_engine& generator)
    {
        gradients.resize(PERLIN_TABLE_SIZE);
        permutations.resize(PERLIN_TABLE_SIZE);

        std::uniform_real_distribution<RealType> distribution(-1.0, 1.0);

        //create our permutations table and stuff
        for (unsigned int i = 0; i < PERLIN_TABLE_SIZE; i++)
        {
            permutations[i] = i;

            gradients[i] = Vector2(distribution(generator), distribution(generator));
            gradients[i].normalize();
        }

        std::shuffle(permutations.begin(), permutations.end(), generator);
    }

public:

    /// Initialize with a random seed
    PerlinGenerator()
    {
        std::default_random_engine generator(std::random_device{}());
        initialize(generator);
    }

    /// Initialize with a given seed
    PerlinGenerator(std::size_t seed)
    {
        std::default_random_engine generator(seed);
        initialize(generator);
    }

    /// evaluate single octave of gradient noise with given location and frequency
    RealType gradientNoise(Vector2 atPos, RealType frequency)
    {
        Vector2 dArrayIndices = atPos * frequency;

        Vector2 baseIndices = Vector2( std::floor(dArrayIndices[0]), std::floor(dArrayIndices[1]));

        Vector2 lerpParams = dArrayIndices - baseIndices;

        Eigen::Vector2i tableIndices =  baseIndices.cast<int>(); //cast to int

        tableIndices[0] &= 0xff;
        tableIndices[1] &= 0xff;

        Eigen::Vector2i incrementedIndices = tableIndices;
        incrementedIndices[0]++; incrementedIndices[1]++;

        incrementedIndices[0] &= 0xff;
        incrementedIndices[1] &= 0xff;

        const int topLeftIndex = (tableIndices[0] + permutations[incrementedIndices[1]]) & 0xff;
        const int topRightIndex = (incrementedIndices[0] + permutations[incrementedIndices[1]]) & 0xff;
        const int bottomLeftIndex = (tableIndices[0] + permutations[tableIndices[1]])&0xff;
        const int bottomRightIndex = (incrementedIndices[0] + permutations[tableIndices[1]])&0xff;

        const Vector2& upperLeftGrad = gradients[topLeftIndex];
        const Vector2& upperRightGrad = gradients[topRightIndex];
        const Vector2& lowerLeftGrad = gradients[bottomLeftIndex];
        const Vector2& lowerRightGrad = gradients[bottomRightIndex];

        Vector2 tempVec4 = lerpParams;
        Vector2 tempVec1 = Vector2(lerpParams[0], -(1.0 - lerpParams[1]));
        Vector2 tempVec2 = Vector2(- (1.0 - lerpParams[0]), -(1.0 - lerpParams[1]));
        Vector2 tempVec3 = Vector2(- (1.0 - lerpParams[0]), lerpParams[1]);

        RealType ulVal = (tempVec1).dot(upperLeftGrad);
        RealType urVal = (tempVec2).dot(upperRightGrad);
        RealType lrVal = (tempVec3).dot(lowerRightGrad);
        RealType llVal = (tempVec4).dot(lowerLeftGrad);

        RealType lerpXTop =  s_interpolate(ulVal, urVal, lerpParams[0]);
        RealType lerpXBottom = s_interpolate(llVal, lrVal, lerpParams[0]);

        RealType result = s_interpolate(lerpXBottom, lerpXTop, lerpParams[1]);

        result = std::min(RealType(1.0), result);
        result = std::max(RealType(-1.0), result);

        return result;
    }

    template<ValidFunctor<RealType> Functor, bool normalize=true>
    RealType fractalSumNoise
    (
        Vector2 atPos,
        int octaves,
        RealType baseFrequency,
        RealType persistence = RealType(.5),
        const Functor& f = [](const RealType& noiseVal) {return noiseVal; }
    )
    {
        RealType frequency = baseFrequency;
        RealType geoAmplitude = 1.0;
        RealType sum = 0.0;

        RealType normFactor = 0.0;

        for (int i = 0; i < octaves; i++, geoAmplitude *= persistence, frequency *= RealType(2.0))
        {
            // assumes that the functor is monotonically increasing so 1.0 from the gradient will map to the largest value in our functor-noise
            normFactor += geoAmplitude * f(RealType(1.0));

            RealType noiseVal = geoAmplitude * f(gradientNoise(atPos, frequency));
            sum += noiseVal;
        }

        if (normalize) sum /= normFactor;

        return sum;
    }

    RealType fractalNoiseAbs(Vector2 atPos, int octaves, RealType baseFrequency, RealType persistence = RealType(.5))
    {
        return fractalSumNoise
        (
            atPos,
            octaves,
            baseFrequency,
            persistence,
            [](const RealType& noiseVal)
            {
                return std::abs<RealType>(noiseVal);
            }
        );
    }

    template<std::uint32_t axis = 0>
    RealType fractalNoiseSin(Vector2 atPos, int octaves, RealType baseFrequency, RealType persistence = RealType(.5))
    {
        static_assert(axis < 2, "fractalNoiseSin() Axis element out of bounds");

        RealType offset = atPos[axis];
        return fractalSumNoise
        (
            atPos,
            octaves,
            baseFrequency,
            persistence,
            [=](const RealType& noiseVal)
            {
                return std::sin<RealType>(offset + fractalNoiseAbs(atPos, octaves, baseFrequency, persistence));
            }
        );
    }
};
