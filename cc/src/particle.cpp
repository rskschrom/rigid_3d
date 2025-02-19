#include "particle.h"
#include "matrix.h"

void Particle::initialize()
{
    // calulate the inertia momentum tensor and its principal axes
    Eigen::Matrix3f matInermPA, eigVecs;
    
    // remove mean from points and calculate inertia moment tensor
    relPoints = centerPoints(points);
    matInerm = inertiaMomentTensor();
    
    std::cout << "inertia moment tensor" << std::endl;
    std::cout << matInerm << std::endl;
    
    // rotate and inertia momentum tensor particle points to principal axes reference frame
    eigVecs = sortedEigVecs(matInerm);
    
    std::cout << "eigenvectors" << std::endl;
    std::cout << eigVecs << std::endl;
    
    matInermPA = eigVecs.transpose() * matInerm * eigVecs;
    
    std::cout << "rotated inertia moment tensor" << std::endl;
    std::cout << matInermPA << std::endl;
    
    relPoints = pointsToBodyFrame(relPoints, eigVecs);
    setMatInerm(matInermPA);
}

float Particle::totalMass()
{
    return pointMass*float(points.size())/3.;
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
