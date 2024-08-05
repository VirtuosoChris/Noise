# Noise

Build status : You should just be able to include this header and use it if you're using the Eigen library.

I just did a bunch of revisions on the 2012 version to make it more type robust and use functor-based post processing of the fractal noise values, but I haven't tested these changes yet since no one is using this code and I did it on a lark.

If you have a problem raise an issue or just use the old header.

You can find-replace for your own vector types if you aren't.  I might make it more generic in the future to work with other vector types with template magic, but it's not a priority.

References:

[Making Noise by Ken Perlin](https://web.archive.org/web/20071011035810/http://noisemachine.com/talk1/)

[Perlin Noise - Wikipedia](https://en.wikipedia.org/wiki/Perlin_noise)
