#bind layer src? val=0
#bind layer !&Cdout

#bind parm cam_rotation float val=20
#bind parm cam_distance float val=10
#bind parm cam_height float val = 1.3

#bind parm defocus_angle float val=0.3
#bind parm focus_dist float val=14

#bind parm samples float val=16
#bind parm max_bounces int val=6

#bind parm num_spheres int val=200
#bind parm spacing float val=1

#import "random.h"

const float PI = 3.14159265;
const int sphere_limit = 5000;

float3 random_in_unit_sphere(float2 p) {
    float3 rand = VEXrandom_2_3(p.x*672.0f, p.y*551.0f);
    
    float phi = 2.0 * PI * rand.x;
    float cosTheta = 2.0 * rand.y - 1.0;
    float u = rand.z;

    float theta = acos(cosTheta);
    float r = pow(u, 1.0f / 3.0f);

    float x = r * sin(theta) * cos(phi);
    float y = r * sin(theta) * sin(phi);
    float z = r * cos(theta);

    return (float3)(x, y, z);
}

float3 random_unit_vector(float2 p) {
    return normalize(random_in_unit_sphere(p));
}

float3 random_in_unit_disk(float2 p) {
    float3 rand = VEXrandom_2_3(p.x*278.0f, p.y*222.0f);
    return (float3)(rand.x*2-1, rand.y*2-1, 0);
}

struct ray {
    float3 origin;
    float3 dir;
};

const int material_lambertian = 0;
const int material_metal = 1;
const int material_dielectric = 2;

struct material {
    int type;
    float3 albedo;
    float metal_fuzz;
    float ior;
};

struct hit_record {
    float3 p;
    float3 normal;
    float t;
    struct material material;
    bool hit;
};

struct sphere {
    float3 center;
    float radius;
    struct material material;
};

struct hit_record hit_sphere(struct sphere sph, struct ray r, struct hit_record rec, bool hit_anything){
    float closest_so_far = rec.t;
    float3 oc = r.origin - sph.center;
    float a = dot(r.dir, r.dir);
    float half_b = dot(oc, r.dir);
    float c = dot(oc, oc) - sph.radius * sph.radius;
    float discriminant = half_b * half_b - a * c;
    
    rec.hit = hit_anything;
    
    if (discriminant < 0.0) {
        return rec;
    }

    float sqrtd = sqrt(discriminant);

    float root = (-half_b - sqrtd) / a;
    if (root < 0.001 || closest_so_far < root) {
        root = (-half_b + sqrtd) / a;
        if (root < 0.001 || closest_so_far < root) {
            return rec;
        }
    }
    
    float3 p = r.origin + r.dir * root;
    
    rec.p = p;
    rec.normal = (p - sph.center) / sph.radius;
    rec.t = root;
    rec.material = sph.material;
    rec.hit = true;
    
    return rec;
}

struct sphere spheres[sphere_limit];

struct hit_record hit(struct ray r, struct hit_record rec, int sphere_count) {
    bool hit = false;
    rec.hit = false;
    
    rec.p = (float3)(0.0f, 0.0f, 0.0f);
    rec.normal = (float3)(0.0f, 0.0f, 0.0f);
    rec.t = 9999.0f;
    rec.material.type = material_lambertian;
    rec.material.albedo = (float3)(0.0f, 0.0f, 1.0f);
    rec.material.metal_fuzz = 0.0;
    rec.material.ior = 0.;
    
    for( int i = 0; i < sphere_count; i++){
        rec = hit_sphere(spheres[i], r, rec, hit);
        hit = rec.hit;
    }
    return rec;
}

bool near_zero(float3 p) {
    float s = 1e-8;
    return p.x < s && p.y < s && p.z < s;
}

float reflectance(float cosine, float ref_idx) {
    float r0 = (1. - ref_idx) / (1. + ref_idx);
    r0 = r0 * r0;
    return r0 + (1. - r0) * pow((1. - cosine), 5.);
}

float3 reflect( float3 v, float3 n){
    return v - 2*dot(v,n)*n;
}

float3 refract(float3 uv, float3 n, float etai_over_etat) {
    float cos_theta = fmin(dot(-uv, n), 1.0f);
    float3 r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    
    float3 r_out_parallel = -sqrt(fabs(1.0f - pow(length(r_out_perp), 2))) * n;
    return r_out_perp + r_out_parallel;
}

// splitting scatter function to scatter_ray and scatter_attenuation
struct ray scatter_ray(struct hit_record rec, struct ray r, float2 seed, float3 attenuation, struct ray scattered){
    struct material m = rec.material;
    
    float3 scatter_direction = normalize(rec.normal + random_unit_vector(seed));
    
    if (m.type == material_lambertian) {
        float3 scatter_direction = normalize(rec.normal + random_unit_vector(seed));

        if (near_zero(scatter_direction)) {
           scatter_direction = rec.normal;
        }
    }
    if (m.type == material_metal) {
        float3 reflected = reflect(r.dir, rec.normal);
        
        struct ray scattered_;
        scattered_.origin = rec.p;
        scattered_.dir = normalize(reflected + m.metal_fuzz * random_in_unit_sphere(seed));
        
        if (dot(scattered_.dir, rec.normal) > 0.) {
            scatter_direction = scattered_.dir;
        }
    }
    if (m.type == material_dielectric){
        bool front_face = dot(r.dir, rec.normal) < 0.0f;
        float3 adjusted_normal = front_face ? rec.normal : -rec.normal;
        float ref = m.ior;
        float refraction_ratio = front_face ? 1.0f/ref : ref;

        float cos_theta = min(dot(-r.dir, adjusted_normal), 1.0f);
        float sin_theta = sqrt(1.0f - cos_theta * cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0f;
        float3 direction;

        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > VEXrandom_2_1(seed.x, seed.y)) {
            scatter_direction = reflect(r.dir, adjusted_normal);
        } else {
            scatter_direction = refract(r.dir, adjusted_normal, refraction_ratio);
        }
    }
    
    scattered.origin = rec.p;
    scattered.dir = scatter_direction;
        
    return scattered;
}

float3 scatter_attenuation(struct hit_record rec, struct ray r, float2 seed, float3 attenuation, struct ray scattered){
    struct material m = rec.material;
    
    if (m.type == material_lambertian) {
        attenuation = m.albedo;
    }
    if (m.type == material_metal) {
        attenuation = m.albedo;
    }
    if (m.type == material_dielectric){
        attenuation = (float3){1.0f, 1.0f, 1.0f};
    }
    return attenuation;
}

float3 ray_color(struct ray r, float2 seed, int max_bounces, int sphere_count) {
    float3 color = (float3){1.0f, 1.0f, 1.0f};
    
    struct hit_record rec;
    int depth;
    
    for (depth = 0; depth < max_bounces; depth++) {
        
        rec = hit(r, rec, sphere_count);
    
        if (rec.hit) {
            
            struct ray scattered;
            float3 attenuation;
            
            attenuation = scatter_attenuation(rec, r, seed * 999.0f + (float)depth, attenuation, scattered);
            scattered = scatter_ray(rec, r, seed * 999.0f + (float)depth, attenuation, scattered);
            
            r = scattered;
            color *= attenuation;
            
        } else {
            float t = 0.5 * (r.dir.y + 1.0);
            color *= mix((float3){1.0f, 1.0f, 1.0f}, (float3){0.5f, 0.7f, 1.0f}, t);
            t = 0.5 * (r.dir.x+r.dir.y + 0.75);

            break;
        }
    }

    if (depth == max_bounces) {
        return (float3){0.0f, 0.0f, 0.0f};
    }
    
    return color;
}
    
float mod(float a, float m){
    return a - m*floor(a/m);
}

struct material ground01 = {0, (float3){0.2f, 0.2f , 0.2f}};
struct material dielectric01 = {2, (float3){0.0f, 0.0f , 0.0f}, 0.0f, 1.5f};
struct material diffuse01 = {0, (float3){0.4f, 0.2f , 0.1f}};
struct material metal01 = {1, (float3){0.7f, 0.6f , 0.5f}, 0.5f};

float4 col;
struct material materials[sphere_limit];
float3 sphere_albedo;

@KERNEL
{   
    //ground
    spheres[0].center = (float3)(0.0f, -3000.0f, 0.0f);
    spheres[0].radius = 3000.0f;
    spheres[0].material = ground01;
    
    float rows = sqrt((float)@num_spheres);
    
    int sphere_count = 2;
    
    //scattered spheres
    for( int i = 5; i < @num_spheres; i++){
        float x_pos = (mod(i, rows)*@spacing) - (@spacing*rows)/2;
        float y_pos = ((i/rows)*@spacing) - (@spacing*rows)/2;
        
        x_pos = x_pos + (float)(VEXrandom_1_1(i+949))*@spacing - (@spacing/2);
        y_pos = y_pos + (float)(VEXrandom_1_1(i+1195))*@spacing - (@spacing/2);
        
        if( length((float3){x_pos, 0.2f, y_pos}-(float3){0, 0.2, 0}) > 1.2f
            && length((float3){x_pos, 0.2f, y_pos}-(float3){-4, 0.2, 0}) > 1.2f
            && length((float3){x_pos, 0.2f, y_pos}-(float3){4, 0.2, 0}) > 1.2f){
            float choose_mat = (float)(VEXrandom_1_1(i+584));
            
            spheres[i].center = (float3)(x_pos, 0.2f, y_pos);
            spheres[i].radius = 0.2f;
            
            if( choose_mat < 0.8){
                sphere_albedo = VEXrandom_2_3(i+983, i+734);
                materials[i].type = 0;
                materials[i].albedo = sphere_albedo*sphere_albedo;
                spheres[i].material = materials[i];
            }
            else if( choose_mat < 0.95){
                sphere_albedo = VEXrandom_2_3(i+983, i+734);
                materials[i].type = 1;
                materials[i].metal_fuzz = 0.2f;
                materials[i].albedo = sphere_albedo*sphere_albedo;
                spheres[i].material = materials[i];
            }
            else{
                materials[i].type = 2;
                materials[i].ior = 1.5f;
                spheres[i].material = materials[i];
            }
            sphere_count++;
        }
    }
    
    spheres[1].center = (float3)(0.0f, 1.0f, 0.0f);
    spheres[1].radius = 1.0f;
    spheres[1].material = dielectric01;
    
    spheres[2].center = (float3)(-4.0f, 1.0f, 0.0f);
    spheres[2].radius = 1.0f;
    spheres[2].material = diffuse01;
    
    spheres[3].center = (float3)(4.0f, 1.0f, 0.0f);
    spheres[3].radius = 1.0f;
    spheres[3].material = metal01;
    
    //camera
    
    float2 pindex = (float2){ (float)@ixy.x, (float)@ixy.y};
    float2 res = (float2){ (float)@res.x, (float)@res.y};

    float3 lookfrom = (float3){0.0f, @cam_height, 0.0f};
    float cam_rotation = (@cam_rotation-45.0)*PI/180.0;
    
    lookfrom.x = (cos(cam_rotation) - sin(cam_rotation))*@cam_distance;
    lookfrom.z = (sin(cam_rotation) + cos(cam_rotation))*@cam_distance;
    
    float3 lookat = (float3){0.0f, 0.5f, 0.0f};
    float3 vup = (float3){0.0f, 1.0f, 0.0f};
    float vfov = 20.0f;
    float aspect_ratio = (float)@res.x / (float)@res.y;
    float focus_dist = @focus_dist;

    float theta = radians(vfov);
    float h = tan(theta / 2.0f);
    float viewport_height = 2.0f * h * focus_dist;
    float viewport_width = aspect_ratio * viewport_height;

    float3 w = normalize(lookfrom - lookat);
    float3 u = normalize(cross(vup, w));
    float3 v = cross(w, u);
    
    float3 viewport_u = viewport_width * u;
    float3 viewport_v = viewport_height * v;
    float3 pixel_delta_u = viewport_u / res.x;
    float3 pixel_delta_v = viewport_v / res.y;
    
    float3 viewport_lower_left = lookfrom - (focus_dist * w) - viewport_u/2 - viewport_v/2;
    float3 pixel00_loc = viewport_lower_left + 0.5f * (pixel_delta_u + pixel_delta_v);
    
    float3 center = lookfrom;
    
    float defocus_radius = focus_dist * tan(radians(@defocus_angle/2.0f));
    float3 defocus_disk_u = u * defocus_radius;
    float3 defocus_disk_v = v * defocus_radius;

    // render
    float3 color = (float3){0.0f, 0.0f, 0.0f};
    for (float s = 0.0f; s < @samples; s++) {
        
        float2 rand = VEXrandom_2_2(pindex.x * 999.0f + s, pindex.y * 729.0f + s);

        struct ray r;
        
        //sample square
        float3 offset = (float3)(rand.x-0.5f, rand.y-0.5f, 0.0f);
        
        float3 pixel_sample = pixel00_loc + ((pindex.x + offset.x) * pixel_delta_u) + ((pindex.y + offset.y) * pixel_delta_v);
    
        //defocus disc sample
        float3 p = random_in_unit_disk(rand);
        float3 dds = center + (p.x * defocus_disk_u) + (p.y * defocus_disk_v);
        
        r.origin = (@defocus_angle <= 0) ? center : dds;
        r.dir = normalize(pixel_sample - r.origin);
        
        color += ray_color(r, rand, @max_bounces, sphere_count+2);
    }
    
    col = (float4)(sqrt(color / @samples), 1.0f);
    @Cdout.set(col);
}