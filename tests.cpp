#include "PerlinNoise.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <iostream>

inline std::uint8_t saturatePixel(double x)
{
    double scaledValue = ((x + 1.0) * 0.5) * 255.0; // map [-1, 1] to [0, 255]
    return static_cast<std::uint8_t>(std::clamp(scaledValue, 0.0, 255.0));
}

inline void writePixel(std::uint8_t* buffer, int x, int y, int width, int height, int channels, double val)
{
    int pixelCoord = (x + y * height);

    for (int i = 0; i < channels; i++)
    {
        buffer[pixelCoord * channels + i] = saturatePixel(val);
    }
}

template <typename f>
void makeNoiseImage(const std::string& path, int width, int height, const f& pix)
{
    PerlinGenerator<double> p;

    std::vector<std::uint8_t> imgBuffer;

    const int channels = 3;

    imgBuffer.resize(width * height * channels);

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            double val = pix(p, Eigen::Vector2d(x, y));
            writePixel(imgBuffer.data(), x, y, width, height, channels, val);
        }
    }

    stbi_write_jpg(path.c_str(), width, height, channels, imgBuffer.data(), 100);
}


int main(void)
{
    const auto& grad = [](const PerlinGenerator<double>& p, const Eigen::Vector2d& coord) { return p.gradientNoise(coord, 1.0 / 32.0); };

    const auto& fract = [](const PerlinGenerator<double>& p, const Eigen::Vector2d& coord)
        {
            return p.fractalSumNoise(coord, 4,  1.0 / 256.0);
        };

    const auto& fract2 = [](const PerlinGenerator<double>& p, const Eigen::Vector2d& coord)
        {
            return p.fractalSumNoise(coord, 4, 1.0 / 256.0, .8, [](const double& noiseVal) {return noiseVal; });
        };

    const auto& fractAbs = [](const PerlinGenerator<double>& p, const Eigen::Vector2d& coord)
        {
            return p.fractalNoiseAbs(coord, 4, 1.0 / 64.0);
        };

    const auto& fractRidged = [](const PerlinGenerator<double>& p, const Eigen::Vector2d& coord)
        {
            return p.fractalNoiseRidged(coord, 4, 1.0 / 24.0, .5, 3.0);
        };

    const auto& fractSin = [](const PerlinGenerator<double>& p, const Eigen::Vector2d& coord)
        {
            // original 'marble texture' = imageData2[w + h * height] = 255 * ((sin( .0075*double(w)    + 8.0*3.14*summation )+1.0)*.5);
            double baseFrequency = 1.0 / 32.0;
            double noiseFrequency = 1.0 / 64.0;
            double amplitude = 4.0 * 3.14;
            double offset = coord[0];
            double x = (offset * (2.0 * 3.14259) * baseFrequency);
            return p.fractalNoiseSin(x, amplitude, coord, 8, noiseFrequency);
        };

    std::cout << "Making test1 (gradient slice)" << std::endl;
    makeNoiseImage("gradient_slice.jpg", 256, 256, grad);
    std::cout << "Making test 2 (fractal noise)" << std::endl;
    makeNoiseImage("fractal_sum.jpg", 512, 512, fract);
    std::cout << "Making test 3 (fractal noise, high persistence)" << std::endl;
    makeNoiseImage("fractal_sum_highp.jpg", 512, 512, fract2);

    std::cout << "Making test 4 (abs noise)" << std::endl;
    makeNoiseImage("abs_noise.jpg", 128, 128, fractAbs);


    std::cout << "Making test 4 (sin noise)" << std::endl;
    makeNoiseImage("sin_noise.jpg", 128, 128, fractSin);


    std::cout << "Making test 5 (ridged multifractal noise)" << std::endl;
    makeNoiseImage("ridged_multifractal.jpg", 128, 128, fractRidged);

}