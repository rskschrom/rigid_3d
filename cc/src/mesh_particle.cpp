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
    
    p111 = 2*(alp(1)*alp(2)*bet(0)+alp(0)*alp(2)*bet(1)+alp(0)*alp(1)*bet(2));
    
    // sum to get integral over triangle
    integral = 2.*area*(p300*m300+p030*m030+p003*m003+
                        p120*m120+p210*m210+
                        p102*m102+p201*m201+
                        p012*m012+p021*m021+p111*m111);
    return integral;
}

      
Eigen::Matrix3f MeshParticle::inertiaMomentTensor()
{

    Eigen::Matrix3f inerm = Eigen::Matrix3f::Zero();
    int nface = faces.rows();
    int iv1, iv2, iv3;
    Eigen::Vector3f v1, v2, v3, v21, v31, vcross;
    Eigen::Vector3f x, y, z;
    float txx, tyy, tzz, txy, txz, tyz;
    
    // loop over faces
    for(int i=0; i < nface; i++)
    {
        // get face vertices
        iv1 = faces(i,0);
        iv2 = faces(i,1);
        iv3 = faces(i,2);
        
        v1 = vertices.row(iv1);
        v2 = vertices.row(iv2);
        v3 = vertices.row(iv3);
        
        // get x,y,z arrays
        x(0) = v1(0);
        x(1) = v2(0);
        x(2) = v3(0);
        
        y(0) = v1(1);
        y(1) = v2(1);
        y(2) = v3(1);
        
        z(0) = v1(2);
        z(1) = v2(2);
        z(2) = v3(2);
                
        // calculate triangle integrals
        txx = triAlp2BetIntegral(x, x, faceAreas(i));
        tyy = triAlp2BetIntegral(y, y, faceAreas(i));
        tzz = triAlp2BetIntegral(z, z, faceAreas(i));
        txy = triAlp2BetIntegral(x, y, faceAreas(i));
        txz = triAlp2BetIntegral(x, z, faceAreas(i));
        tyz = triAlp2BetIntegral(y, z, faceAreas(i));
                
        // add to inertia moment tensors
        inerm(0,0) = inerm(0,0)+tyy*faceNorms(i,1)+tzz*faceNorms(i,2);
        inerm(1,1) = inerm(1,1)+txx*faceNorms(i,0)+tzz*faceNorms(i,2);
        inerm(2,2) = inerm(2,2)+txx*faceNorms(i,0)+tyy*faceNorms(i,1);
        
        inerm(0,1) = inerm(0,1)+txy*faceNorms(i,0);
        inerm(0,2) = inerm(0,2)+txz*faceNorms(i,0);
        inerm(1,2) = inerm(1,2)+tyz*faceNorms(i,1);
        
    }
    
    // scale tensor elements and get symmetric elements
    inerm(0,0) = density/3.*inerm(0,0);
    inerm(1,1) = density/3.*inerm(1,1);
    inerm(2,2) = density/3.*inerm(2,2);
    
    inerm(0,1) = density/2.*inerm(0,1);
    inerm(0,2) = density/2.*inerm(0,2);
    inerm(1,2) = density/2.*inerm(1,2);
        
    inerm(1,0) = inerm(0,1);
    inerm(2,0) = inerm(0,2);
    inerm(2,1) = inerm(1,2);
    
    return inerm;
}
/*
void MeshParticle::initialize()
{
    return;
}
*/