#include "HODLR_Node.hpp"
#include <Rcpp.h>

HODLR_Node::HODLR_Node(int level_number, int node_number, int local_number, 
                       int n_start, int n_size, double tolerance
                      ) 
{
    // Storing the passed parameters as the attributes of the created object:
    this->level_number = level_number;
    this->node_number  = node_number;
    this->n_start      = n_start;
    this->n_size       = n_size;
    // Start location and size for the child on the left:
    this->c_start[0]   = n_start;
    this->c_size[0]    = 0.5 * n_size;
    // Start location and size for the child on the right:
    this->c_start[1]   = n_start + c_size[0];
    this->c_size[1]    = n_size  - c_size[0];
    this->tolerance    = tolerance;
}

void HODLR_Node::assembleLeafNode(HODLR_Matrix* A) 
{
    // At the leaf level we are just going to be building the matrix
    // directly since it's a full rank block:
    K = A->getMatrix(n_start, n_start, n_size, n_size);
}

void HODLR_Node::matmatProductLeaf(Mat x, Mat& b) 
{
    b.block(n_start, 0, n_size, x.cols()) += K * x.block(n_start, 0, n_size, x.cols());
}

void HODLR_Node::assembleNonLeafNode(HODLR_Matrix* A, bool is_sym) 
{
    if(is_sym == true)
    {
        A->rookPiv(c_start[0], c_start[1], c_size[0], c_size[1], tolerance, U[0], V[1], rank[0]);
        V[0]    = U[0];
        U[1]    = V[1];
        rank[1] = rank[0];
    }

    else
    {
        A->rookPiv(c_start[0], c_start[1], c_size[0], c_size[1], tolerance, U[0], V[1], rank[0]);
        A->rookPiv(c_start[1], c_start[0], c_size[1], c_size[0], tolerance, U[1], V[0], rank[1]);
    }
}

void HODLR_Node::matmatProductNonLeaf(Mat x, Mat& b) 
{
    b.block(c_start[0], 0, c_size[0], x.cols()) += 
    (U[0] * (V[1].transpose() * x.block(c_start[1], 0, c_size[1], x.cols())));

    b.block(c_start[1], 0, c_size[1], x.cols()) += 
    (U[1] * (V[0].transpose() * x.block(c_start[0], 0, c_size[0], x.cols())));
}

void HODLR_Node::printNodeDetails()
{
	Rcpp::Rcout << "Level Number       :" << level_number << std::endl;
	Rcpp::Rcout << "Node Number        :" << node_number << std::endl;
	Rcpp::Rcout << "Start of Node      :" << n_start <<  std::endl;
	Rcpp::Rcout << "Size of Node       :" << n_size <<  std::endl;
	Rcpp::Rcout << "Tolerance          :" << tolerance << std::endl;

    for(int i = 0; i < 2; i++)
    {
        if(i == 0)
            Rcpp::Rcout << "Left Child:" << std::endl;
        else
            Rcpp::Rcout << "Right Child:" << std::endl;

	Rcpp::Rcout << "Start of Child Node:" << c_start[i] <<  std::endl;
	Rcpp::Rcout << "Size of Child Node :" << c_size[i]  <<  std::endl;
    }

    Rcpp::Rcout << "Shape of U[0]      :" << U[0].rows() << ", " << U[0].cols() << std::endl;
    Rcpp::Rcout << "Shape of U[1]      :" << U[1].rows() << ", " << U[1].cols() << std::endl;
    Rcpp::Rcout << "Shape of V[0]      :" << V[0].rows() << ", " << V[0].cols() << std::endl;
    Rcpp::Rcout << "Shape of V[1]      :" << V[1].rows() << ", " << V[1].cols() << std::endl;
    Rcpp::Rcout << "Shape of K         :" << K.rows() << ", " << K.cols() << std::endl;
}
