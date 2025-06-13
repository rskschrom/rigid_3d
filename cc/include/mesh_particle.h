#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include "io.h"

#ifndef MESH_PARTICLE_H
#define MESH_PARTICLE_H

// barycentric integral constants
const float m300 = 1./20.;
const float m030 = 1./20.;
const float m003 = 1./20.;
const float m120 = 1./60.;
const float m210 = 1./60.;
const float m012 = 1./60.;
const float m021 = 1./60.;
const float m102 = 1./60.;
const float m201 = 1./60.;
const float m111 = 1./120.;

/// Class for triangulated mesh particle object
class MeshParticle
{

    private:
        Eigen::Matrix3f matInerm;

    public:
    
        /// initial properties
        Eigen::MatrixX3f vertices;
        Eigen::MatrixX3i faces;
        float density;
        
        // computed properties of the mesh
        Eigen::VectorXf faceAreas;
        Eigen::MatrixX3f faceNorms;

        /*!
         * Default constructor.
         *
         * \param vertices the vertices of the mesh.
         * \param faces the indices of each face corresponding to a triangular patch.
         * \param density the density of the particle.
         */    

        MeshParticle(Eigen::MatrixX3f vertices, Eigen::MatrixX3i faces, float density):
                     vertices(vertices),
                     faces(faces),
                     density(density) { initialize(); }

        /*!
         * Initializing by reorienting the particle to the principal axes of its inertia momentum tensor.
         *
         */
         void initialize();
         
         /*!
         * Calculate the area and normal vector of each face.
         *
         */
         void calculateFaceAreasNorms();

        /*!
         * Calculate the total particle mass.
         *
         * \return the total particle mass.
         */
        float totalMass();
        
        /*!
         * Calculate the particle inertia momentum tensor.
         *
         * \return the particle inertia momentum tensor.
         */
        Eigen::Matrix3f inertiaMomentTensor();
        
        float triAlp2BetIntegral(Eigen::Vector3f alpPoints, Eigen::Vector3f betPoints, float area);
        
        void setMatInerm(Eigen::Matrix3f imat){ matInerm = imat; }
        
        Eigen::Matrix3f getMatInerm() { return matInerm; }
        
        Eigen::MatrixX3f getVertices() { return vertices; }
        
        Eigen::VectorXf getFaceAreas() { return faceAreas; }
        
        Eigen::MatrixX3f getFaceNorms() { return faceNorms; }
                
        /*!
         * Calculate the particle center of mass.
         *
         * \return the particle center of mass.
         */
        std::vector<float> centerOfMass();
        
        /*!
         * Translate the particle.
         *
         * \param tr the translation vector.
         */
        void translate(std::vector<float> tr);
        
        /*!
         * Rotate the particle.
         *
         * \param quat the rotation quaternion.
         */
        void rotate(Eigen::Matrix3f rmat);
        
        

};

#endif
