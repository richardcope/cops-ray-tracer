# Ray Tracing in One Cops Node

A simple ray tracer, written to function as an OpenCL snippet within the COPS context of SideFX Houdini.\
Created from the concepts described in the indispensable Ray Tracing in One Weekend book:\
https://raytracing.github.io/books/RayTracingInOneWeekend.html

![alt text](https://github.com/richardcope/cops-ray-tracer/blob/main/images/capture001.JPG "rendered scene in cops")

## To Use

1) Drop down an OpenCl node within a Cops context
2) Paste cops_ray_tracer.c into the kernel code text box
3) Set signature>outputs1 to Cdout and type to RGB
4) Click icon above kernel code to add parameters to interface\

There is also an example hip file in this repository.\
\
\
https://www.youtube.com/watch?v=JjPsg5_ESIU


