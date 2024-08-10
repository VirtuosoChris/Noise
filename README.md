# Noise

Build status : You should just be able to include this header and use it if you're using the Eigen library.  
You can find-replace for your own vector types if you aren't using Eigen.  I might make it more generic in the future to work with other vector types with template magic, but it's not a priority.
There's a test case that won't work on its own unless you add eigen and stb_image_write.  There's a folder in the repo with the output of that test program.

I made a bunch of revisions on the 2012 version to make it more type robust and use functor-based post processing of the fractal noise values.  I got some notes from Steve132 on that revision and I think the end result is pretty nice.

I want to extend it to have a shader generator, N-dimensions, and multicore.

# References:

[Making Noise by Ken Perlin](https://web.archive.org/web/20071011035810/http://noisemachine.com/talk1/)

[Perlin Noise - Wikipedia](https://en.wikipedia.org/wiki/Perlin_noise)

# Sample Images
![ridged_multifractal](https://github.com/user-attachments/assets/d27c59b6-3de5-4bc6-9f60-81b0ea39fd81)
![test5_8pi](https://github.com/user-attachments/assets/8e4a53fa-44a8-4a06-8016-182d3a0a52cb)
![test5](https://github.com/user-attachments/assets/213864b5-117d-4a7c-864b-1022e700a0b7)
![test4](https://github.com/user-attachments/assets/3a7a69ea-26fd-4671-a0f0-96a464985760)
