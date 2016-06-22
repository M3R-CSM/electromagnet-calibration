#ifndef EIGENTOYAML_H
#define EIGENTOYAML_H
/**
    This header file creates template specializations to work with the libyaml-cpp library for Eigen matricies and std::vector<double> vectors.

**/


#include<eigen3/Eigen/Core>
#include<yaml-cpp/yaml.h>
#include<vector>

#include <sstream>
using std::stringstream;
#include <iostream>
using namespace std;

namespace YAML {
template<>
struct convert<Eigen::VectorXd> {
    static Node encode(const Eigen::VectorXd& rhs) {
        Node node;
        for( unsigned int i=0; i<rhs.size(); i++)
        {
            node.push_back( rhs(i) );
        }
        return node;
    }

    static bool decode(const Node& node, Eigen::VectorXd& rhs) {

        if( node.size() == 0 && node.IsSequence()  ) {
            // Empty Sequence
            rhs.resize(0,1);
            return true;

        }else if( node.IsSequence() && node[0].size() == 1 ) {
            rhs.resize(node.size(),1);

            for( unsigned int i=0; i<node.size(); i++ )
                rhs[i] =  node[i][0].as<double>();

            return true;

        } else if( node.IsSequence() && node[0].size() == 0 ) {
            rhs.resize(node.size(),1);

            for( unsigned int i=0; i<node.size(); i++ )
                rhs[i] =  node[i].as<double>() ;

            return true;
        } else  if( node.IsScalar() )
        {
            stringstream tmp;
            tmp << node;
            string nodeStr = tmp.str();

            if( (nodeStr.size()==1  && nodeStr.c_str()[0] == '~') )
                rhs.resize(0,1);
            else
            {
                rhs.resize(1,1);
                rhs[0] = node.as<double>();
            }

            return true;
        }

        /*if(!node.IsSequence() || ( node.size() > 0 && node[0].size() != 1) )
        {
            return false;
        }

        if( node.size() == 0 )
        {
            rhs.resize(0,0);
            return true;
        }

        rhs.resize(node.size(),1);

        for( unsigned int i=0; i<node.size(); i++ )
            rhs[i] = node[i][0].as<double>();

        return true;*/
    }
};

template<>
struct convert<Eigen::Vector3d> {
    static Node encode(const Eigen::Vector3d& rhs) {
        Node node;
        for( unsigned int i=0; i<3; i++)
        {
            node.push_back(rhs(i));
        }
        return node;
    }

    static bool decode(const Node& node, Eigen::Vector3d& rhs) {
        if(!node.IsSequence() || node.size() != 3 ) {
            return false;
        }

        if( node[0].size() == 1 )
            for( unsigned int i=0; i<node.size(); i++ )
                rhs[i] = node[i][0].as<double>();
        else if( node[0].size() == 0 )
            for( unsigned int i=0; i<node.size(); i++ )
                rhs[i] = node[i].as<double>();
        else
            return false;


        return true;
    }
};

template<>
struct convert<Eigen::MatrixXd> {
    static Node encode(const Eigen::MatrixXd& rhs) {
        Node node;
        for( unsigned int i=0; i<rhs.rows(); i++)
        {
            /*Node tmpNode;
            for( unsigned int j=0; j<rhs.cols(); j++ )
                tmpNode.push_back( rhs(i,j) );

            node.push_back( tmpNode );
            */
            for( unsigned int j=0; j<rhs.cols(); j++ )
                node[i][j] = rhs(i,j);
        }
        return node;
    }

   static bool decode(const Node& node, Eigen::MatrixXd& rhs) {
        if( !node.IsSequence() ) {
            return false;
        }

	if( node.size() == 0 ) 
	{
            // Empty Sequence
            rhs.resize(0,0);
            return true;
        }

        rhs.resize( node.size(), node[0].size() );

        for( unsigned int i=0; i<node.size(); i++ )
            for( unsigned int j=0; j<node[0].size(); j++ )
                rhs(i,j) = node[i][j].as<double>();

        return true;
    }
};

template<>
struct convert<Eigen::Matrix3d> {
    static Node encode(const Eigen::Matrix3d& rhs) {
        Node node;
        for( unsigned int i=0; i<rhs.rows(); i++)
        {
            //Node tmpNode;
            for( unsigned int j=0; j<rhs.cols(); j++ )
                //tmpNode.push_back( rhs(i,j) );
                node[i][j] = rhs(i,j);

            //node.push_back( tmpNode );
        }
        return node;
    }

   static bool decode(const Node& node, Eigen::Matrix3d& rhs) {
        if( !node.IsSequence() )
        {
            return false;
        }

        if( node.size() != 3 || node[0].size() != 3 )
	{
            // Not a 3x3 matrix!
            return false;
        }

        rhs.resize( node.size(), node[0].size() );

        for( unsigned int i=0; i<node.size(); i++ )
            for( unsigned int j=0; j<node[0].size(); j++ )
                rhs(i,j) = node[i][j].as<double>();

        return true;
    }
};



template<>
struct convert< std::vector<double> > {
    static Node encode(const std::vector<double>& rhs) {
        Node node;
        for( unsigned int i=0; i<rhs.size(); i++)
        {
            node.push_back((rhs[i]));
        }
        return node;
    }

    static bool decode(const Node& node, std::vector<double>& rhs) {

        if( node.size() == 0) // && node.IsSequence()
        {
            // Empty Sequence
            rhs.resize(0);
            return true;

        }else if( node.IsSequence() && node[0].size() == 1 ) {
            for( unsigned int i=0; i<node.size(); i++ )
                rhs.push_back( node[i][0].as<double>() );

            return true;

        } else if( node.IsSequence() && node[0].size() == 0 ) {
            for( unsigned int i=0; i<node.size(); i++ )
                rhs.push_back( node[i].as<double>() );

            return true;
        } else  if( node.IsScalar() )
        {
            stringstream tmp;
            tmp << node;
            string nodeStr = tmp.str();

            if( (nodeStr.size()==1  && nodeStr.c_str()[0] == '~') )
                rhs.resize(0);
            else
                rhs.push_back(node.as<double>());

            return true;
        }

        // Otherwise
        return false;
    }
};

}

#endif // EIGENTOYAML_H
