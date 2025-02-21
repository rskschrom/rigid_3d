#include "mesh_particle.h"
#include "matrix.h"

float MeshParticle::totalMass()
{
    std::cout << vertices.rows() << " " << vertices.cols() << std::endl;
    return float(vertices.rows());
}

void MeshParticle::calculateFaceAreas()
{
    int nface = faces.rows();
    int iv1, iv2, iv3;
    Eigen::Vector3f v1, v2, v3, v21, v31, vcross;
    float area;
    
    // initialize face area vector
    faceAreas.resize(nface);
    
    for(int i=0; i < nface; i++)
    {
        // get face vertices
        iv1 = faces(i,0);
        iv2 = faces(i,1);
        iv3 = faces(i,2);
        
        v1 = vertices.row(iv1);
        v2 = vertices.row(iv2);
        v3 = vertices.row(iv3);
        
        // calculate (positive) area with cross product
        v21 = v2-v1;
        v31 = v3-v1;
        vcross = v21.cross(v31);
        faceAreas(i) = 0.5*sqrt(vcross(0)*vcross(0)+vcross(1)*vcross(1)+vcross(2)*vcross(2));
        
        std::cout << v1 << " " << faceAreas(i) << std::endl;
        
    }
    
    std::cout << faces.rows() << " " << faces.cols() << std::endl;
    return;
}

/*
void MeshParticle::initialize()
{
    return;
}
        
Eigen::Matrix3f Particle::inertiaMomentTensor()
{

    Eigen::Matrix3f inerm = Eigen::Matrix3f::Zero();
    int npar = relPoints.size()/3;
    
    // loop over tensor axes
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
      
            // tensor off diagonal
            if (i!=j){
                for (int k = 0; k < npar; k++){
                    inerm(i,j) += -relPoints[3*k+i]*relPoints[3*k+j];
                }
            }
      
            // tensor diagonal
            else{
                for (int k = 0; k < npar; k++){
                    inerm(i,j) += relPoints[3*k]*relPoints[3*k]+
                                  relPoints[3*k+1]*relPoints[3*k+1]+
                                  relPoints[3*k+2]*relPoints[3*k+2]-
                                  relPoints[3*k+i]*relPoints[3*k+i];
                }
            }
            
            inerm(i,j) = pointMass*inerm(i,j);
        }
    }
    return inerm;
}
*/