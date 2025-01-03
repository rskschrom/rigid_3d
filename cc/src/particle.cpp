#include "particle.h"
#include "matrix.h"

float Particle::totalMass()
{
    return pointMass*float(relPoints.size())/3.;
}

std::vector<float> Particle::transformPoints()
{
    int npar = relPoints.size()/3;
    float a, b, c;
    std::vector<float> prot(3);
    std::vector<float> absPos(npar*3);

    // get rotation coefficients
    a = orient[0]*orient[0]-
        orient[1]*orient[1]-
        orient[2]*orient[2]-
        orient[3]*orient[3];
      
    c = -2.*orient[0];
  
    // loop over points
    for (int i = 0; i < npar; i++){
      b = 2.*(relPoints[3*i]*orient[1]+
              relPoints[3*i+1]*orient[2]+
              relPoints[3*i+2]*orient[3]);
            
      // prot = a*p+b*v+c*pxv
      prot[0] = a*relPoints[3*i]+b*orient[1]+
                c*(relPoints[3*i+1]*orient[3]-relPoints[3*i+2]*orient[2]);
              
      prot[1] = a*relPoints[3*i+1]+b*orient[2]+
                c*(relPoints[3*i+2]*orient[1]-relPoints[3*i]*orient[3]);
              
      prot[2] = a*relPoints[3*i+2]+b*orient[3]+
                c*(relPoints[3*i]*orient[2]-relPoints[3*i+1]*orient[1]);
    
      absPos[3*i] = prot[0]+comPos[0];
      absPos[3*i+1] = prot[1]+comPos[1];
      absPos[3*i+2] = prot[2]+comPos[2];
    }
  
    return absPos;
  
}
        
std::vector<float> Particle::inertiaMomentTensor()
{

    std::vector<float> inerm(9);
    int npar = relPoints.size()/3;
    
    // loop over tensor axes
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
      
            // tensor off diagonal
            if (i!=j){
                for (int k = 0; k < npar; k++){
                    inerm[3*i+j] += -relPoints[3*k+i]*relPoints[3*k+j];
                }
            }
      
            // tensor diagonal
            else{
                for (int k = 0; k < npar; k++){
                    inerm[3*i+j] += relPoints[3*k]*relPoints[3*k]+
                                    relPoints[3*k+1]*relPoints[3*k+1]+
                                    relPoints[3*k+2]*relPoints[3*k+2]-
                                    relPoints[3*k+i]*relPoints[3*k+i];
                }
            }
            
            inerm[3*i+j] = pointMass*inerm[3*i+j];
        }
    }
    return inerm;
}

void Particle::initialize()
{
    // get inertia tensor and inverse
    std::vector<float> inerm = inertiaMomentTensor();
    Eigen::Matrix3f matInerm = fvecMat(inerm);
    Eigen::Matrix3f matIInerm = matrixInv(matInerm);
    
    // get principal axes of inertia momentum tensor
    Eigen::Matrix3f matInermPA, matIInermPA;
    Eigen::EigenSolver<Eigen::Matrix3f> es(matInerm);
    Eigen::Matrix3f eigVecs = es.eigenvectors().real();
    Eigen::Matrix3f permute = Eigen::Matrix3f::Zero();
    
    // sort eigenvectors by decreasing eigenvalue
    float a = es.eigenvalues().real()[0];
    float b = es.eigenvalues().real()[1];
    float c = es.eigenvalues().real()[2];
    float tmpVal;
    int tmpInd;
    std::vector<int> indSort = {0, 1, 2};
    std::vector<float> eigVals = {a, b, c};
    
    //std::cout << eigVals[0] << "\t" << eigVals[1] << "\t" << eigVals[2] << std::endl;
    
    if (eigVals[0]>eigVals[1]){
        tmpVal = eigVals[0];
        tmpInd = indSort[0];
        eigVals[0] = eigVals[1];
        eigVals[1] = tmpVal;
        indSort[0] = indSort[1];
        indSort[1] = tmpInd;
    }
    if (eigVals[0]>eigVals[2]){
        tmpVal = eigVals[0];
        tmpInd = indSort[0];
        eigVals[0] = eigVals[2];
        eigVals[2] = tmpVal;
        indSort[0] = indSort[2];
        indSort[2] = tmpInd;
    }
    if (eigVals[1]>eigVals[2]){
        tmpVal = eigVals[1];
        tmpInd = indSort[1];
        eigVals[1] = eigVals[2];
        eigVals[2] = tmpVal;
        indSort[1] = indSort[2];
        indSort[2] = tmpInd;
    }
    
    //std::cout << indSort[0] << "\t" << indSort[1] << "\t" << indSort[2] << std::endl;
    
    // permute matrix
    for (int i = 0; i<3; i++){
        permute(i,indSort[i]) = 1.;
    }
    
    eigVecs = eigVecs * permute;
    
    matInermPA = eigVecs.transpose() * matInerm * eigVecs;
    matIInermPA = eigVecs.transpose() * matIInerm * eigVecs;
    setMatInerm(matInermPA);
    setMatIInerm(matIInermPA);
    
}

void Particle::write()
{
    writeVector(comPos, "position.txt");
    writeVector(orient, "orientation.txt");
}

