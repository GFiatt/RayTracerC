# RayTracer in C

## 📄 Description
**RayTracer in C** renders 3D scenes made of **spheres** and **polygons** (n-gons) using the Phong model (**ambient, diffuse, specular**), **distance attenuation**, **hard shadows**, and **mirror reflections** (recursive).  
The program **reads a scene from a text file** and displays the image in an OpenGL window; it also saves `output.jpg`.

**Specs**
- Resolution: **1008 × 567** (`WIDTH`, `HEIGHT`)
- Up to **6** point lights
- **4× supersampling**
- Outputs: `output.avs` (raw) and `output.jpg` (via **libjpeg**)

---

## 🚀 How to Run

### Linux (Debian/Ubuntu)
```bash
sudo apt-get update
sudo apt-get install build-essential freeglut3-dev libjpeg-dev
make
./PR04 poolTable.txt
```

### Windows (MSYS2 / MinGW)
1) Install MSYS2 and open **MSYS2 MinGW x64**.  
2) Dependencies:
```bash
pacman -S --needed mingw-w64-x86_64-toolchain mingw-w64-x86_64-freeglut mingw-w64-x86_64-libjpeg-turbo
```
3) Build and run:
```bash
gcc -O2 raytracer.c -o PR04.exe -lfreeglut -lopengl32 -lglu32 -ljpeg -lm
./PR04.exe assets/poolTable.txt
```

> You can pass any scene file: `./PR04 myScene.txt`.

---

## 🛠️ Features
- **Primitives:** spheres and **polygons** (N≥3, coplanar).
- **Lighting:** Phong (Ka, Kd, Ks, Kn) + per-light `Ip` + **attenuation** `1/(a + b·d²)`.
- **Shadows:** occlusion checks toward each light.
- **Reflections:** blended with **Kr** (limited recursion depth).
- **Anti-aliasing:** 4 samples per pixel.
- **I/O:** scene from `.txt`, OpenGL/GLUT display, and **JPEG** export.

---

## ✍️ Scene File Format

Plain-text file; blocks in **any order**: `LIGHT`, `SPHERE`, `POLYGON`.  
Values are `long double`. Comments are not supported.

### LIGHT
```
LIGHT
x= <x>
y= <y>
z= <z>
Ip= <intensity>      # typically 0..1
```

### SPHERE
```
SPHERE
xc= <center_x>
yc= <center_y>
zc= <center_z>
r= <radius>
R= <0..1>            # RGB color
G= <0..1>
B= <0..1>
Ka= <ambient>
Kd= <diffuse>
Ks= <specular>
Kn= <shininess>      # ~10..200
Kr= <reflection>     # 0..1
```

### POLYGON
```
POLYGON
borders= <N>=3..n
x= <x0>, y= <y0>, z= <z0>
x= <x1>, y= <y1>, z= <z1>
...
x= <xN-1>, y= <yN-1>, z= <zN-1>
R= <0..1>
G= <0..1>
B= <0..1>
Ka= <ambient>
Kd= <diffuse>
Ks= <specular>
Kn= <shininess>
Kr= <reflection>
```

**Notes**
- Polygons must be **coplanar**; the normal is computed from two edges and normalized.
- Point-in-polygon test is performed on the **dominant-axis projection** of the polygon’s plane (ray-crossing).
- Global ambient light is fixed in code: `Ia = 0.8`.

---

## 🧪 Example Input

Run:
```
./PR04 poolTable.txt
```

---

## 📁 Project Structure
```
.
├─ poolTable.txt
├─ PR04.c
├─ README.md
└─ Makefile
```

---

## 👨‍💻 Contributors
- Jose Andrés Ramírez  
- **Gabriel Fiatt** (maintainer)

---

## 📎 License
MIT License — see `LICENSE`.
