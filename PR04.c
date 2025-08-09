#include <GL/glut.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <jpeglib.h>
#include <string.h>

#define WIDTH 1008
#define HEIGHT 567
#define EPSILON 0.00005

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Vector {
    long double x;
    long double y;
    long double z;
};

struct RGB {
    long double red;
    long double green;
    long double blue;
};

struct Sphere {
    long double Xc;
    long double Yc;
    long double Zc;
    long double r;
    struct RGB color;
    long double Ka;
    long double Kd;
    long double Ks;
    long double Kn;
    long double Kr;
};

struct Polygon {
    struct Vector *vertices;
    int num_vertices;
    struct Vector normal;
    struct RGB color;
    long double Ka;
    long double Kd;
    long double Ks;
    long double Kn;
    long double Kr;
};

struct Light {
    struct Vector position;
    long double Ip;
};

struct Intersection {
    long double t;
    int is_sphere; 
    struct Sphere *sphere;
    struct Polygon *polygon;
};

struct RGB Background = {0.45, 0.45, 0.45};
struct Vector eye_Vector = {0.0, 40.0, -200.0};
int NUM_SPHERES;
struct Sphere *spheres;

int NUM_POLYGONS;
struct Polygon *polygons;

int NUM_LIGHTS = 6;
struct Light lights[6];

struct RGB FrameBuffer[WIDTH][HEIGHT];

long double Ia = 0.8;

void leer_escena_desde_archivo(const char* nombre_archivo) {
    FILE *archivo = fopen(nombre_archivo, "r");
    if (!archivo) {
        printf("No se pudo abrir el archivo de entrada\n");
        exit(1);
    }

    NUM_SPHERES = 0;
    NUM_POLYGONS = 0;
    NUM_LIGHTS = 0;
    char tipo[10];
    char buffer[256];
    while (fscanf(archivo, "%s", tipo) != EOF) {
        if (strcmp(tipo, "LIGHT") == 0) {
            if (NUM_LIGHTS < 6) {
                fscanf(archivo, " x= %Lf", &lights[NUM_LIGHTS].position.x);
                fscanf(archivo, " y= %Lf", &lights[NUM_LIGHTS].position.y);
                fscanf(archivo, " z= %Lf", &lights[NUM_LIGHTS].position.z);
                fscanf(archivo, " Ip= %Lf", &lights[NUM_LIGHTS].Ip);
                NUM_LIGHTS++;
            }
        } else if (strcmp(tipo, "POLYGON") == 0) {
            polygons = realloc(polygons, (NUM_POLYGONS + 1) * sizeof(struct Polygon));
            int num_vertices;
            fscanf(archivo, " borders= %d", &num_vertices);
            polygons[NUM_POLYGONS].num_vertices = num_vertices;
            polygons[NUM_POLYGONS].vertices = malloc(num_vertices * sizeof(struct Vector));
            for (int j = 0; j < num_vertices; j++) {
                fscanf(archivo, " x= %Lf, y= %Lf, z= %Lf",
                       &polygons[NUM_POLYGONS].vertices[j].x,
                       &polygons[NUM_POLYGONS].vertices[j].y,
                       &polygons[NUM_POLYGONS].vertices[j].z);
            }
            fscanf(archivo, " R= %Lf", &polygons[NUM_POLYGONS].color.red);
            fscanf(archivo, " G= %Lf", &polygons[NUM_POLYGONS].color.green);
            fscanf(archivo, " B= %Lf", &polygons[NUM_POLYGONS].color.blue);
            fscanf(archivo, " Ka= %Lf", &polygons[NUM_POLYGONS].Ka);
            fscanf(archivo, " Kd= %Lf", &polygons[NUM_POLYGONS].Kd);
            fscanf(archivo, " Ks= %Lf", &polygons[NUM_POLYGONS].Ks);
            fscanf(archivo, " Kn= %Lf", &polygons[NUM_POLYGONS].Kn);
            fscanf(archivo, " Kr= %Lf", &polygons[NUM_POLYGONS].Kr);
            
            struct Vector edge1 = {
                polygons[NUM_POLYGONS].vertices[2].x - polygons[NUM_POLYGONS].vertices[0].x,
                polygons[NUM_POLYGONS].vertices[2].y - polygons[NUM_POLYGONS].vertices[0].y,
                polygons[NUM_POLYGONS].vertices[2].z - polygons[NUM_POLYGONS].vertices[0].z
            };
            struct Vector edge2 = {
                polygons[NUM_POLYGONS].vertices[1].x - polygons[NUM_POLYGONS].vertices[0].x,
                polygons[NUM_POLYGONS].vertices[1].y - polygons[NUM_POLYGONS].vertices[0].y,
                polygons[NUM_POLYGONS].vertices[1].z - polygons[NUM_POLYGONS].vertices[0].z
            };
            polygons[NUM_POLYGONS].normal.x = edge1.y * edge2.z - edge1.z * edge2.y;
            polygons[NUM_POLYGONS].normal.y = edge1.z * edge2.x - edge1.x * edge2.z;
            polygons[NUM_POLYGONS].normal.z = edge1.x * edge2.y - edge1.y * edge2.x;
            long double length = sqrt(polygons[NUM_POLYGONS].normal.x * polygons[NUM_POLYGONS].normal.x +
                                      polygons[NUM_POLYGONS].normal.y * polygons[NUM_POLYGONS].normal.y +
                                      polygons[NUM_POLYGONS].normal.z * polygons[NUM_POLYGONS].normal.z);
            polygons[NUM_POLYGONS].normal.x /= length;
            polygons[NUM_POLYGONS].normal.y /= length;
            polygons[NUM_POLYGONS].normal.z /= length;
            NUM_POLYGONS++;
        } else if (strcmp(tipo, "SPHERE") == 0) {
            spheres = realloc(spheres, (NUM_SPHERES + 1) * sizeof(struct Sphere));
            fscanf(archivo, " xc= %Lf", &spheres[NUM_SPHERES].Xc);
            fscanf(archivo, " yc= %Lf", &spheres[NUM_SPHERES].Yc);
            fscanf(archivo, " zc= %Lf", &spheres[NUM_SPHERES].Zc);
            fscanf(archivo, " r= %Lf", &spheres[NUM_SPHERES].r);
            fscanf(archivo, " R= %Lf", &spheres[NUM_SPHERES].color.red);
            fscanf(archivo, " G= %Lf", &spheres[NUM_SPHERES].color.green);
            fscanf(archivo, " B= %Lf", &spheres[NUM_SPHERES].color.blue);
            fscanf(archivo, " Ka= %Lf", &spheres[NUM_SPHERES].Ka);
            fscanf(archivo, " Kd= %Lf", &spheres[NUM_SPHERES].Kd);
            fscanf(archivo, " Ks= %Lf", &spheres[NUM_SPHERES].Ks);
            fscanf(archivo, " Kn= %Lf", &spheres[NUM_SPHERES].Kn);
            fscanf(archivo, " Kr= %Lf", &spheres[NUM_SPHERES].Kr);
            NUM_SPHERES++;
        }
    }

    fclose(archivo);
}

int calculate_intersection(struct Vector ancla, struct Vector dir, struct Sphere *sphere, long double *t) {
    long double Xd = dir.x;
    long double Yd = dir.y;
    long double Zd = dir.z;
    long double xe = ancla.x;
    long double ye = ancla.y;
    long double ze = ancla.z;
    long double Xc = sphere->Xc;
    long double Yc = sphere->Yc;
    long double Zc = sphere->Zc;
    long double r = sphere->r;

    long double a = Xd*Xd + Yd*Yd + Zd*Zd;
    long double b = 2.0 * (Xd*(xe - Xc) + Yd*(ye - Yc) + Zd*(ze - Zc));
    long double c = (xe - Xc)*(xe - Xc) + (ye - Yc)*(ye - Yc) + (ze - Zc)*(ze - Zc) - r*r;

    long double discriminant = b*b - 4.0*a*c;

    if (discriminant < 0) {
        return 0;
    } else {
        long double sqrt_discriminant = sqrt(discriminant);
        long double t0 = (-b - sqrt_discriminant) / (2.0*a);
        long double t1 = (-b + sqrt_discriminant) / (2.0*a);

        if (t0 > 0.001 && t1 > 0.001) {
            *t = fmin(t0, t1);
            return 1;
        } else if (t0 > 0.001) {
            *t = t0;
            return 1;
        } else if (t1 > 0.001) {
            *t = t1;
            return 1;
        } else {
            return 0;
        }
    }
}

long double producto_punto(struct Vector N, struct Vector L) {
    return N.x * L.x + N.y * L.y + N.z * L.z;
}

struct Intersection first_intersection(struct Vector ancla, struct Vector dir) {
    struct Intersection closest_intersection;
    closest_intersection.t = INFINITY;
    closest_intersection.is_sphere = -1;
    closest_intersection.sphere = NULL;
    closest_intersection.polygon = NULL;

    for (int i = 0; i < NUM_SPHERES; i++) {
        long double t;
        if (calculate_intersection(ancla, dir, &spheres[i], &t)) {
            if ((t < closest_intersection.t) && (t > EPSILON)) {
                closest_intersection.t = t;
                closest_intersection.is_sphere = 1;
                closest_intersection.sphere = &spheres[i];
            }
        }
    }

    for (int i = 0; i < NUM_POLYGONS; i++) {
        long double t;
        struct Vector P;
        struct Vector N = polygons[i].normal;
        long double denom = producto_punto(N, dir);
        if (fabs(denom) > EPSILON) {
            long double d = producto_punto(N, polygons[i].vertices[0]);
            t = (d - producto_punto(N, ancla)) / denom;
            if (t > EPSILON && t < closest_intersection.t) {
                P.x = ancla.x + dir.x * t;
                P.y = ancla.y + dir.y * t;
                P.z = ancla.z + dir.z * t;
                int crossings = 0;
                int proj_axis = 0;
                struct Vector N_abs = {fabs(N.x), fabs(N.y), fabs(N.z)};
                if (N_abs.x > N_abs.y && N_abs.x > N_abs.z) {
                    proj_axis = 0;
                } else if (N_abs.y > N_abs.z) {
                    proj_axis = 1;
                } else {
                    proj_axis = 2;
                }
                for (int j = 0; j < polygons[i].num_vertices; j++) {
                    struct Vector v1 = polygons[i].vertices[j];
                    struct Vector v2 = polygons[i].vertices[(j + 1) % polygons[i].num_vertices];
                    long double u0, v0, u1, v1_p, uP, vP;
                    switch (proj_axis) {
                        case 0:
                            u0 = v1.y; v0 = v1.z;
                            u1 = v2.y; v1_p = v2.z;
                            uP = P.y; vP = P.z;
                            break;
                        case 1:
                            u0 = v1.x; v0 = v1.z;
                            u1 = v2.x; v1_p = v2.z;
                            uP = P.x; vP = P.z;
                            break;
                        case 2:
                            u0 = v1.x; v0 = v1.y;
                            u1 = v2.x; v1_p = v2.y;
                            uP = P.x; vP = P.y;
                            break;
                    }
                    if (((v0 > vP) != (v1_p > vP)) && fabs(v1_p - v0) > EPSILON) {
                        long double interseccion = (u1 - u0) * (vP - v0) / (v1_p - v0) + u0;
                        if (uP < interseccion) {
                            crossings++;
                        }
                    }
                }
                if (crossings % 2 == 1) {
                    closest_intersection.t = t;
                    closest_intersection.is_sphere = 0;
                    closest_intersection.polygon = &polygons[i];
                }
            }
        }
    }

    return closest_intersection;
}

struct RGB de_que_color(struct Vector ancla, struct Vector dir, int depth) {
    struct RGB color;
    if (depth > 3) {
        color = Background;
        return color;
    }

    struct Intersection intersection = first_intersection(ancla, dir);
    if (intersection.is_sphere == -1) {
        color = Background;
    } else {
        long double tmin = intersection.t;
        struct Vector P = {ancla.x + dir.x * tmin, ancla.y + dir.y * tmin, ancla.z + dir.z * tmin};
        struct Vector N;
        struct RGB obj_color;
        long double Ka, Kd, Ks, Kn, Kr;

        if (intersection.is_sphere) {
            N.x = (P.x - intersection.sphere->Xc) / intersection.sphere->r;
            N.y = (P.y - intersection.sphere->Yc) / intersection.sphere->r;
            N.z = (P.z - intersection.sphere->Zc) / intersection.sphere->r;
            obj_color = intersection.sphere->color;
            Ka = intersection.sphere->Ka;
            Kd = intersection.sphere->Kd;
            Ks = intersection.sphere->Ks;
            Kn = intersection.sphere->Kn;
            Kr = intersection.sphere->Kr;
        } else {
            N = intersection.polygon->normal;
            obj_color = intersection.polygon->color;
            Ka = intersection.polygon->Ka;
            Kd = intersection.polygon->Kd;
            Ks = intersection.polygon->Ks;
            Kn = intersection.polygon->Kn;
            Kr = intersection.polygon->Kr;
        }

        long double N_len = sqrt(N.x * N.x + N.y * N.y + N.z * N.z);
        N.x /= N_len;
        N.y /= N_len;
        N.z /= N_len;

        struct Vector V = { -dir.x, -dir.y, -dir.z };
        long double I = Ia * Ka;
        long double E = 0.0;

        for (int k = 0; k < NUM_LIGHTS; k++) {
            struct Vector L = {lights[k].position.x - P.x, lights[k].position.y - P.y, lights[k].position.z - P.z};
            long double L_len = sqrt(L.x * L.x + L.y * L.y + L.z * L.z);
            L.x /= L_len;
            L.y /= L_len;
            L.z /= L_len;

            long double dot = producto_punto(N, L);
            if (dot > EPSILON) {
                struct Vector R = {2 * N.x * dot - L.x, 2 * N.y * dot - L.y, 2 * N.z * dot - L.z};
                long double dot2 = producto_punto(R, V);

                struct Vector P_offset = {P.x + N.x * EPSILON, P.y + N.y * EPSILON, P.z + N.z * EPSILON};
                long double Fatt = 1.0 / (0.1 + 0.00002 * L_len * L_len);

                struct Intersection obstaculo = first_intersection(P_offset, L);

                int mismo_objeto = 0;
                if (obstaculo.is_sphere == intersection.is_sphere) {
                    if (obstaculo.is_sphere && obstaculo.sphere == intersection.sphere) {
                        mismo_objeto = 1;
                    } else if (!obstaculo.is_sphere && obstaculo.polygon == intersection.polygon) {
                        mismo_objeto = 1;
                    }
                }

                if (obstaculo.is_sphere == -1 || obstaculo.t > L_len - EPSILON || mismo_objeto) {
                    I += lights[k].Ip * Kd * dot * Fatt;
                    if (dot2 > EPSILON) {
                        E += pow(dot2, Kn) * Ks * lights[k].Ip * Fatt;
                    }
                }
            }
        }

        I = fmin(I, 1.0);
        E = fmin(E, 1.0);

        struct RGB local_color;
        local_color.red = (obj_color.red * I) + E * (1.0 - obj_color.red);
        local_color.green = (obj_color.green * I) + E * (1.0 - obj_color.green);
        local_color.blue = (obj_color.blue * I) + E * (1.0 - obj_color.blue);

        if (Kr > EPSILON) {
            struct Vector R = {
                dir.x - 2 * producto_punto(dir, N) * N.x,
                dir.y - 2 * producto_punto(dir, N) * N.y,
                dir.z - 2 * producto_punto(dir, N) * N.z
            };

            struct Vector P_offset = {P.x + R.x * EPSILON, P.y + R.y * EPSILON, P.z + R.z * EPSILON};
            struct RGB reflected_color = de_que_color(P_offset, R, depth + 1);

            color.red = (1 - Kr) * local_color.red + Kr * reflected_color.red;
            color.green = (1 - Kr) * local_color.green + Kr * reflected_color.green;
            color.blue = (1 - Kr) * local_color.blue + Kr * reflected_color.blue;
        } else {
            color = local_color;
        }
    }

    return color;
}

void guardar_avs(const char *nombre_archivo) {
    FILE *archivo = fopen(nombre_archivo, "wb");
    if (!archivo) {
        printf("No se pudo abrir el archivo para escribir\n");
        return;
    }

    unsigned int width = WIDTH;
    unsigned int height = HEIGHT;
    fwrite(&width, sizeof(unsigned int), 1, archivo);
    fwrite(&height, sizeof(unsigned int), 1, archivo);

    for (int j = 0; j < HEIGHT; j++) {
        for (int i = 0; i < WIDTH; i++) {
            unsigned char red = (unsigned char)(FrameBuffer[i][j].red *255);
            unsigned char green = (unsigned char)(FrameBuffer[i][j].green *255);
            unsigned char blue = (unsigned char)(FrameBuffer[i][j].blue *255);
            fwrite(&red, sizeof(unsigned char), 1, archivo);
            fwrite(&green, sizeof(unsigned char), 1, archivo);
            fwrite(&blue, sizeof(unsigned char), 1, archivo);
        }
    }

    fclose(archivo);
}

void leer_avs(const char *nombre_archivo, unsigned char **pixel_data, int *width, int *height) {
    FILE *archivo = fopen(nombre_archivo, "rb");
    if (!archivo) {
        printf("No se pudo abrir el archivo AVS\n");
        return;
    }
    fread(width, sizeof(unsigned int), 1, archivo);
    fread(height, sizeof(unsigned int), 1, archivo);

    *pixel_data = (unsigned char *)malloc(3 * (*width) * (*height));

    fread(*pixel_data, sizeof(unsigned char), 3 * (*width) * (*height), archivo);

    fclose(archivo);
}

void guardar_jpeg(const char *nombre_jpeg, unsigned char *pixel_data, int width, int height) {
    FILE *archivo = fopen(nombre_jpeg, "wb");
    if (!archivo) {
        printf("No se pudo abrir el archivo JPEG\n");
        return;
    }

    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    jpeg_stdio_dest(&cinfo, archivo);

    cinfo.image_width = width;
    cinfo.image_height = height;
    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_RGB;

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, 100, TRUE);
    jpeg_start_compress(&cinfo, TRUE);

    JSAMPROW row_pointer;
    int row_stride = width * 3;

    while (cinfo.next_scanline < cinfo.image_height) {
        row_pointer = &pixel_data[cinfo.next_scanline * row_stride];
        jpeg_write_scanlines(&cinfo, &row_pointer, 1);
    }

    jpeg_finish_compress(&cinfo);
    fclose(archivo);
    jpeg_destroy_compress(&cinfo);
}

void rayTracing() {
    int samples = 4;
    long double fov = M_PI / 3.0;
    long double aspect_ratio = (long double)WIDTH / (long double)HEIGHT;
    long double image_plane_distance = 1.0;

    long double image_plane_height = 2.0 * tan(fov / 2.0) * image_plane_distance;
    long double image_plane_width = image_plane_height * aspect_ratio;

    struct Vector image_plane_center = {
        eye_Vector.x,
        eye_Vector.y,
        eye_Vector.z + image_plane_distance
    };

    struct Vector horizontal = {image_plane_width, 0, 0};
    struct Vector vertical = {0, image_plane_height, 0};

    struct Vector lower_left_corner = {
        image_plane_center.x - horizontal.x / 2.0 - vertical.x / 2.0,
        image_plane_center.y - horizontal.y / 2.0 - vertical.y / 2.0,
        image_plane_center.z - horizontal.z / 2.0 - vertical.z / 2.0
    };

    for (int j = 0; j < HEIGHT; j++) {
        for (int i = 0; i < WIDTH; i++) {
            struct RGB color = {0.0, 0.0, 0.0};
            for (int s = 0; s < samples; s++) {
                long double u = ((long double)i + ((long double)rand() / RAND_MAX)) / (long double)(WIDTH - 1);
                long double v = ((long double)(HEIGHT - j - 1) + ((long double)rand() / RAND_MAX)) / (long double)(HEIGHT - 1);

                struct Vector dir = {
                    lower_left_corner.x + u * horizontal.x + v * vertical.x - eye_Vector.x,
                    lower_left_corner.y + u * horizontal.y + v * vertical.y - eye_Vector.y,
                    lower_left_corner.z + u * horizontal.z + v * vertical.z - eye_Vector.z
                };

                long double L = sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
                dir.x /= L;
                dir.y /= L;
                dir.z /= L;

                struct RGB sample_color = de_que_color(eye_Vector, dir, 0);
                color.red += sample_color.red;
                color.green += sample_color.green;
                color.blue += sample_color.blue;
            }

            color.red /= samples;
            color.green /= samples;
            color.blue /= samples;

            FrameBuffer[i][j] = color;
        }
    }

    guardar_avs("output.avs");
    unsigned char *pixel_data;
    int width1, height1;
    leer_avs("output.avs", &pixel_data, &width1, &height1);
    guardar_jpeg("output.jpg", pixel_data, width1, height1);
    free(pixel_data);

    unsigned char pixels[HEIGHT][WIDTH][3];

    for (int i = 0; i < WIDTH; i++) {
        for (int j = 0; j < HEIGHT; j++) {
            pixels[HEIGHT - j - 1][i][0] = (unsigned char)(FrameBuffer[i][j].red * 255);
            pixels[HEIGHT - j - 1][i][1] = (unsigned char)(FrameBuffer[i][j].green * 255);
            pixels[HEIGHT - j - 1][i][2] = (unsigned char)(FrameBuffer[i][j].blue * 255);
        }
    }

    glClear(GL_COLOR_BUFFER_BIT);
    glDrawPixels(WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, pixels);
    glutSwapBuffers();
}

void init(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(WIDTH, HEIGHT);
    glutCreateWindow("Ray Tracing");

    glClearColor(0.45, 0.45, 0.45, 1.0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, WIDTH, 0, HEIGHT);

    glutDisplayFunc(rayTracing);

    glutMainLoop();
}

int main(int argc, char** argv) {
    const char* archivo_entrada = (argc > 1) ? argv[1] : "input.txt";
    leer_escena_desde_archivo(archivo_entrada);
    init(argc, argv);

    free(spheres);
    for (int i = 0; i < NUM_POLYGONS; i++) {
        free(polygons[i].vertices);
    }
    free(polygons);
    return 0;
}

/*
    Wikipedia contributors. (2024, enero 1). Angle of view (photography). Wikipedia. https://en.wikipedia.org/wiki/Angle_of_view_(photography)
*/