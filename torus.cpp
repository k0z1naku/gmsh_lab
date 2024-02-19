#include <set>
#include <gmsh.h>
#include<math.h>
#include<iostream>
#include<vector>

int main(int argc, char **argv)
{
    gmsh::initialize();
    gmsh::model::add("tor");

    const double PI = 3.1415926;
    double R = 10;
    double r1 = 5;
    double r2 = 3;
    double lc = 0.2;
    int N = 50;
    int N_DOTS = N * N;

    //Add all Points
    for (int i = 0; i < N; i++) {
        for (int j = 1; j <= N; j++) {
            double omega = i * 2 * PI / N;
            double phi = j * 2 * PI / N;
            gmsh::model::geo::addPoint((R + r1 * cos(phi))*sin(omega), (R + r1 * cos(phi))*cos(omega), r1 * sin(phi), lc,  i * N + j);
            gmsh::model::geo::addPoint((R + r2 * cos(phi))*sin(omega), (R + r2 * cos(phi))*cos(omega), r2 * sin(phi), lc, i * N + j + N_DOTS);
        }
    }

    //Create points and surface for points of radius r1
    for (int i=0; i<N; i++) {
        for(int j=1; j<=N; j++) {
            gmsh::model::geo::addLine(i * N + j, ((i + 1) % N) * N + j, i * N + j);
        }
    }

    for (int i=0; i<N; i++) {
        for(int j=1; j<=N; j++) {
            gmsh::model::geo::addLine(j + N * i, j % N + 1 + N * i, i * N + j + N_DOTS);
        }
    }

    for(int i=0; i<N; i++) {
        for(int j=1; j<=N; j++) {
            gmsh::model::geo::addCurveLoop({j + i*N, j + N_DOTS + ((i+1)%N)*N, -((j%N+1)+i*N), -(j+N_DOTS + i*N)}, j+i*N);
            gmsh::model::geo::addPlaneSurface({j+i*N}, j+i*N);
        }
    }

    //Create points and surface for points of radius r2
    for (int i=0; i<N; i++) {
        for(int j=1; j<=N; j++) {
            gmsh::model::geo::addLine(i * N + j + N_DOTS, ((i + 1) % N) * N + j + N_DOTS, i * N + j + 2*N_DOTS);
        }
    }

    for (int i=0; i<N; i++) {
        for(int j=1; j<=N; j++) {
            gmsh::model::geo::addLine(j + N * i + N_DOTS, j % N + 1 + N * i + N_DOTS, (i) * N + j + 3* N_DOTS);
        }
    }

    for(int i=0; i<N; i++) {
        for(int j=1; j<=N; j++) {
            gmsh::model::geo::addCurveLoop({j + i*N + 2*N_DOTS, j + 3*N_DOTS + ((i+1)%N)*N, -((j%N+1)+i*N + 2* N_DOTS), -(j+3* N_DOTS + i*N)}, j+i*N + 2*N_DOTS);
            gmsh::model::geo::addPlaneSurface({j+i*N + 2 * N_DOTS}, j+i*N + 2 * N_DOTS);
        }
    }

    //Thicker wall
    for (int i = 0; i < N; i++) {
        std::vector<int> curveTagsWide;
        std::vector<int> curveTagsTight;
        for (int j = 0; j < N; j++) {
            curveTagsWide.push_back(j+1 + i * N + N_DOTS);
            curveTagsTight.push_back(j+1 + i * N + 3*N_DOTS);
        }
        gmsh::model::geo::addCurveLoop(curveTagsWide, (i+1)*N + 3 * N_DOTS);
        gmsh::model::geo::addCurveLoop(curveTagsTight, (i+1)*N + 4*N_DOTS);
        gmsh::model::geo::addPlaneSurface({(i+1)*N +3 *  N_DOTS, -((i+1)*N + 4 * N_DOTS)}, (i+1)*N + N_DOTS);
    }

    //mesh wall
    for(int i=0; i<N; i++) {
        std::vector<int> curveSurfaceTags;
        for(int j=1; j<=N; j++) {
            curveSurfaceTags.push_back(j+i*N);
            curveSurfaceTags.push_back(j+i*N + 2*N_DOTS);
        }
        curveSurfaceTags.push_back((i+1)* N + N_DOTS);
        curveSurfaceTags.push_back(((i+2)%(N+1) + (int) (i+2)/(N+1) )* N + N_DOTS);
        gmsh::model::geo::addSurfaceLoop(curveSurfaceTags, (i+1)*N);
        gmsh::model::geo::addVolume({(i+1)*N}, (i+1)*N);
    }

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(3);
    gmsh::write("mesh/tor.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();

    return 0;
}