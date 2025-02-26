#include "mesh_particle.h"
#include "matrix.h"

float MeshParticle::totalMass()
{
    int nface = faces.rows();
    int iv1, iv2, iv3;
    Eigen::Vector3f centroid, v1, v2, v3, v21, v31;
    float mass = 0.;
    
    for(int i=0; i < nface; i++)
    {
        // get face vertices
        iv1 = faces(i,0);
        iv2 = faces(i,1);
        iv3 = faces(i,2);
        
        v1 = vertices.row(iv1);
        v2 = vertices.row(iv2);
        v3 = vertices.row(iv3);
        
        // calculate centroid
        centroid = (v1+v2+v3)/3.;
        mass = mass+faceAreas(i)/3.*centroid.dot(faceNorms.row(i));
        
    }
    
    return mass*density;
}

void MeshParticle::calculateFaceAreasNorms()
{
    int nface = faces.rows();
    int iv1, iv2, iv3;
    Eigen::Vector3f v1, v2, v3, v21, v31, vcross;
    float area;
    
    // initialize face area vector
    faceAreas.resize(nface);
    faceNorms.resize(nface,3);
    
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
        
        // get normal vector
        //v21mag = sqrt(v21(0)*v21(0)+v21(0)*v21(0));
        faceNorms(i,0) = vcross(0)/(2.*faceAreas(i));
        faceNorms(i,1) = vcross(1)/(2.*faceAreas(i));
        faceNorms(i,2) = vcross(2)/(2.*faceAreas(i));
        
    }
    
    return;
}

float MeshParticle::triAlp2BetIntegral(Eigen::Vector3f alp, Eigen::Vector3f bet, float area)
{
    float integral, p300, p030, p003,
        p120, p210, p102, p201, p012, p021, p111;
    
    // calculate polynomial values
    p300 = alp(0)*alp(0)*bet(0);
    p030 = alp(1)*alp(1)*bet(1);
    p003 = alp(2)*alp(2)*bet(2);
    
    p120 = alp(1)*alp(1)*bet(0)+2.*alp(0)*alp(1)*bet(1);
    p210 = alp(0)*alp(0)*bet(1)+2.*alp(1)*alp(0)*bet(0);
    p102 = alp(2)*alp(2)*bet(0)+2.*alp(0)*alp(2)*bet(2);
    p201 = alp(0)*alp(0)*bet(2)+2.*alp(2)*alp(0)*bet(0);
    p012 = alp(2)*alp(2)*bet(1)+2.*alp(1)*alp(2)*bet(2);
    p021 = alp(1)*alp(1)*bet(2)+2.*alp(2)*alp(1)*bet(1);
    
    p111 = alp(1)*alp(2)*bet(0)+alp(0)*alp(2)*bet(1)+alp(0)*alp(1)*bet(2);
    
    // sum to get integral over triangle
    integral = 2.*area*(p300*m300+p030*m030+p003*m003+
                        p120*m120+p210*m210+
                        p102*m102+p201*m201+
                        p012*m012+p021*m021+p111*m111);
    return integral;
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