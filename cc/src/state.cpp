#include "state.h"
#include "matrix.h"

float State::rotationalPotentialEnergy(float g, float rhofGrad, float rhob){
    Eigen::Matrix3f matInerm = getMatInerm();
    Eigen::Matrix3f matInermWorld, matRot;
    float ixx, iyy, izz;
    float rotPE;
    
    // get inertia mometum tensor in world frame
    matRot = quatToMatrix(orient);
    matInermWorld = matRot * (matInerm * matRot.transpose());
    
    ixx = matInermWorld(0,0);
    iyy = matInermWorld(1,1);
    izz = matInermWorld(2,2);
    
    rotPE = -g*rhofGrad/(4.*rhob)*(ixx+iyy-izz);
    
    return rotPE;
}

float State::rotationalKineticEnergy(){
    Eigen::Matrix3f matInerm = getMatInerm();
    float rotKE;
        
    rotKE = 0.5*(matInerm(0,0)*omega[0]*omega[0]+
                 matInerm(1,1)*omega[1]*omega[1]+
                 matInerm(2,2)*omega[2]*omega[2]);
        
    return rotKE;
}

void State::initialize()
{
    // get inertia tensor and inverse
    Eigen::Matrix3f matInerm = getMatInerm();
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

void State::write()
{
    writeVector(comPos, "position.txt");
    writeVector(orient, "orientation.txt");
}
