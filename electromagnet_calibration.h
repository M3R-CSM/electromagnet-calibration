#ifndef ELECTROMAGNET_CALIBRATION_H
#define ELECTROMAGNET_CALIBRATION_H

/*
 ElectromagnetCalibration is a calss that returns the mangetic field and its various derivatives
    based on a spherical harmonic calibration. The constructor takes a YAML formatted calibration file. An example file is provided below.

    Constructors & Initialization:
        ElectromagnetCalibration( fileName ) : loads a calibration file in YAML format.

        ElectromagnetCalibration( systemName, workSpace, coilList, dc_field_offset ) : Constructs a calibration based off of definition provided.

        loadCalibration( fileName ) : loads a  new calibration file in YAML format.  Returns true if sucessful.

    Field Information Query:
        fieldAtPoint( current, position ) :  Returns the field at the specified point given the currents packed into a vector

        gradientAtPoint( current position) : Returns the 3x3 gradient matrix given the currents packed into a vetctor

        fieldAndGradientAtPoint(current, position) : Returns a 8x1 vector of field stacked on a 5 element vector packing
                of the unique elements in the gradient matrix ([dBx/dx, dBx/dy, dBx/dz, dBy/dy, dBy/dz]).

        offsetFieldAndGradientAtPoint( position) : Returns a 8x1 vector of field stacked on a 5 element vector packing
                of the unique elements in the gradient matrix ([dBx/dx, dBx/dy, dBx/dz, dBy/dy, dBy/dz]) that is due to
                the offset only.

        fieldCurrentJacobian( position ) : Returns a 3xN matrix of the local current to field relationship at a point

        gradientCurrentJacobian( position) : Returns a 5xN matrix of the local current to gradient relationship at a point

        fieldAndGradientCurrentJacobian( position ) : Returns an 8xN matrix of the local current to field stacked ontop of gradient

        gradientPositionJacobian( current, position ) : Returns a 5x3 matrix stating how the gradient tensor of the magnetic field changes with position

        fullMagneticState( fieldAtPoint_out, fieldGradientPositionJacobian_out, fieldGradientCurrentJacobian_out, current_in, position_in ) : uses pass by refence to
                 provide key control parametes about the magetic state at a position.  It is more efficient to call this funciton than to call the above sequentially.

        fullMagneticState( currnet, position ) :  Returns a magnetic state struct that contains the field information at the point requested


    Gradient Manipulation:
        remapGradientVector( gradVector ) : Takes a 5x1 element gradient vector and remaps it to the 3x3 symetiric and zero-trace field gradient matrix

        remapGradientMatrix( gradMatrix) : Takes a 3x3 symetiric and zero-trace field gradient matrix and maps it to a 5x1 vector

        packForceMatrix( moment) :  Takes a 3x1 magnetic dipole moment and packs it into a 3x5 matrix that when multiplied by the 5x1 packing of the gradient matrix yields the force on the dipole.



    Model Parameter Information:
        getNumberOfCoils() : Returns the number of coils (currents) for the systems

        getNumberOfSources( coilNumber ) : Returns the number of ScalorPotentials used to model the field for the coil specified

        getNumberOfCoeffients( coilNumber, potentialNumber ) : Returns the number of coefficients used to model the particular source

        getName() : Returns the Name of the system from the calibration file

        pointInWorkspace( position ) : Returns true if the point is in the calibrated workspace

        useOffset() : Returns if the offset is enabled

        useOffset( true_to_use ) : Enables the offset if true is passed; Disables if false is passed.

        useOffset( newOffset ) : Sets a new ScalorPotential offset and enables it to be used.

        hasOffset() : returns true if an offset has been provided

        getWorkSpace() : returns the workspace structure that reflects the calibrated region

        setWorkSpace( workspace ) : specifices the workspace structure that reflects the calibrated region

    Calibration:
        calibrate(calibrationName, dataList, printProgressIfTrue, printStatisticsIfTrue,
                  calibration_constraints, minimumSourceToCenterDistance, maximumSourceToCenterDistance,
                  converganceTolerance, maxIterations, numberOfConvergedIterations ) : Specifies the data and parameters
                  for calibration and calibrates the potential.

        writeCalibration( fileName ) : writes a YAML formated calibration to file

        printStats(  dataList ) : prints to screen (cout) the calibration statistics



    Libraries:
         Eigen matrix algebra library, version 3 or higher: http://eigen.tuxfamily.org
         YAML-CPP parser and emitter in C++ library: https://github.com/jbeder/yaml-cpp.git

    Files:
        This class requires the following files to compile
            electromagnet_calibration.h
            electromagnet_calibration.cpp
            EigenToYAML.h

    Example calibration file layout:
        System_Name: 'Test System'
        Workspace_Dimensions:
        - [-0.025, 0.025]
        - [-0.025, 0.025]
        - [-0.025, 0.025]
        Coil_List: [Coil_1, Coil_2, Offset]
        Coil_1:
          Source_List: [Src_1, Src_2]
          Src_1:
            A_Coeff: []
            B_Coeff:
            - [1]
            - [-0.5]
            Source_Position:
            - [0.22416258792214033]
            - [-0.07556653205541619]
            - [-0.22109796762044504]
            Source_Direction:
            - [-0.7025088251776764]
            - [0.14854822954743013]
            - [0.6959991249740192]
          Src_2:
            A_Coeff: []
            B_Coeff:
            - [0.25]
            - [-0.25]
            Source_Position:
            - [0.3168094119654824]
            - [0.21831564272785464]
            - [-0.0776279052373666]
            Source_Direction:
            - [-0.96373336339233]
            - [-0.25663768967502687]
            - [0.07317870247555225]
        Coil_2:
          Source_List: [Src_1]
          Src_1:
            A_Coeff: []
            B_Coeff:
            - [1]
            Source_Position:
            - [0.15865180505201498]
            - [-0.022607524623308707]
            - [-0.34998355269929454]
            Source_Direction:
            - [-0.6483538763510639]
            - [0.07075237800794039]
            - [0.7580444261800644]
        Offset:
          Source_List: [Src_1]
          Src_1:
            A_Coeff: []
            B_Coeff:
            - [.01]
            Source_Position:
            - [0.15865180505201498]
            - [-0.022607524623308707]
            - [-0.34998355269929454]
            Source_Direction:
            - [-0.6483538763510639]
            - [0.07075237800794039]
            - [0.7580444261800644]

*/

#include <Eigen/Dense>

#include <assert.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
using std::ofstream;
using std::stringstream;

#include "EigenToYAML.h"
#include "scalorPotential.h"

#include <yaml-cpp/yaml.h>


typedef Eigen::Matrix<double,8,1> Vector8d;
typedef Eigen::Matrix<double,5,1> Vector5d;


/**
 * @brief The MagneticMeasurementData struct provides the format calibration data must be supplied for the calibration
 */
struct MagneticMeasurement
{
    Eigen::Vector3d Field; /**< The measured Field in Tesla. A 3x1 Vector**/
    Eigen::Vector3d Position;/**< The position of the measurement in meters. A 3x1 Vector **/
    Eigen::VectorXd AppliedCurrentVector; /**< The applied current vector in Amps. A Nx1 Vector **/
    MagneticMeasurement(); /**< Initializes all parameters to zero. **/

    /**
     *
     *  \brief calibrationDataPoint
     *  \param Field the field value in Tesla
     *  \param Pos the position of the measurement in meters
     *  \param CurrentVec the applied current vector in Amps
     */
    MagneticMeasurement( const Eigen::Vector3d& Field, const Eigen::Vector3d& Pos, const Eigen::VectorXd& CurrentVec );

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/**
 * @brief The MagneticState struct contains informtion necessary to quantify the field at a position for control.
 */
struct MagneticState
{
    /** @brief Field  The magnetic field at the position in Tesla. This is a 3x1 Matrix. **/
    Eigen::Vector3d Field;

    /** @brief Gradient  The Magnetic gradient at the position in Tesla/Meter.  This is a 3x3 Matrix. **/
    Eigen::Matrix3d Gradient;

    /**
     * @brief GradientPositionJacobian The rate of change of the 5 unique gradient terms with respect to position in Tesla/Meter^2.
     *        This is a 5x3 Matrix, to get the 3x3x3 Gradeint derivative tensor, the function
     *        ElectromagnetCalibration::remapGradientVector can be used on each column to get the gradient matrix
     *        change in x,y and z in that order.
     **/
    Eigen::Matrix<double,5,3> GradientPositionJacobian;

    /**
     * @brief FieldGradientActuationMatrix  The field and gradient (packed into a 5x1 vector) current jacobian.  This is a 8xN matrix,
     *        where N is the number of current sources. The first 3 rows are for field the last five rows are for gradient
     *        ElectromagnetCalibration::remapGradientVector.
     **/
    Eigen::MatrixXd FieldGradientActuationMatrix;

    MagneticState()
    {
        Field.setZero(3);
        Gradient.setZero(3,3);
        GradientPositionJacobian.setZero(5,3);
        FieldGradientActuationMatrix.setZero(0,0);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/**
 *  @brief a class describing the calibration of an arbitrary magnetic system
 *
 *  This class assumes the system can be described by spherical harmonic scalor
 *  potentials at multiple locations with weightings proportional to the currents
 *  appied.
 **/
class ElectromagnetCalibration
{
public:
    /*
     * @brief The MagneticWorkSpace struct defines a rectangular workspace that all of the calibration points lied within to identify if a point is within the region calibrated.
     */
    struct MagneticWorkSpace{
        double xMin; /**< Minimum X Position in Meters */
        double xMax; /**< Maximum X position in Meters */
        double yMin; /**< Minimum Y Position in Meters */
        double yMax; /**< Maximum Y Position in Meters */
        double zMin; /**< Minimum Z Position in Meters */
        double zMax; /**< Maximum Z Position in Meters */
        MagneticWorkSpace();
        MagneticWorkSpace(double size);
        MagneticWorkSpace(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
    };

    /**
     * @brief constructor that builds a calibrated electromagnet system based on a prexisting yaml encoded calibration file.
     *
     * @param calibrationFileName is the file location for the calibration
     **/
    ElectromagnetCalibration( std::string calibrationFileName );

    /**
     *  \brief constructor that builds a calibrated electromanget system based on a list of scalor potentials and an optional offset.
     *  \param workspace_ The workspace for which this calibration applies.
     *  \param coilList The list of coils and their respective sources.
     *  \param dc_field_offset The dc offset field, if any.
     **/
    ElectromagnetCalibration( std::string systemName,  const MagneticWorkSpace& workSpace_, const std::vector<ScalorPotential>& coilList, const ScalorPotential& dc_field_offset = ScalorPotential() );

    // The following functions return the field and/or gradient given a current vector and a position.
    //   If no position is provided, it is assumed to be at the workspace origin. For the combined vector,
    //   the gradient matrix has been repacked into a 5 element vetor form (because it is symetric and has zero trace).
    //   The order of the gradient terms is: [dBx/dx, dBx/dy, dBx/dz, dBy/dy, dBy/dz].

    /**
     * @brief returns the field at a point
     *  @param currentVector is an orderd list of currents in each coil
     *  @param position is the position in the workspace the field is desired
     *
     *  This function does not check to see if the point is actually in the calibrated workspace.
     **/
    Eigen::Vector3d fieldAtPoint( const Eigen::VectorXd& currentVector, const Eigen::Vector3d& position = Eigen::Vector3d::Zero() ) const;

    /**
     * @brief returns the 3x3 symetric gradient matrix at a desired location
     *  @param currentVector is an orderd list of currents in each coil
     *  @param position is the position in the workspace the field is desired
     *
     *  This function does not check to see if the point is actually in the calibrated workspace.
     **/
    Eigen::Matrix3d gradientAtPoint( const Eigen::VectorXd& currentVector, const Eigen::Vector3d& position = Eigen::Vector3d::Zero() ) const;

    /**
     * @brief returns the 8x1 stacked field over gradient vector
     *  @param currentVector is an orderd list of currents in each coil
     *  @param position is the position in the workspace the field is desired
     *
     *  The gradient matrix has been repacked, since it is symetric an has zero trace, into a five element vector.
     *    The element order is: [dBx/dx, dBx/dy, dBx/dz, dBy/dy, dBy/dz]^T
     *  This function does not check to see if the point is actually in the calibrated workspace.
     **/
    Vector8d fieldAndGradientAtPoint( const Eigen::VectorXd& currentVector, const Eigen::Vector3d& position = Eigen::Vector3d::Zero() ) const;

    /**
     * @brief returns the 8x1 stacked field over gradient vector due to the zero-current field offset
     *  @param position is the position in the workspace the field is desired
     *
     *  The gradient matrix has been repacked, since it is symetric an has zero trace, into a five element vector.
     *    The element order is: [dBx/dx, dBx/dy, dBx/dz, dBy/dy, dBy/dz]^T
     *  This function does not check to see if the point is actually in the calibrated workspace.
     **/
    Vector8d offsetFieldAndGradientAtPoint( const Eigen::Vector3d& position = Eigen::Vector3d::Zero() ) const;

    // The following functions return the Current Jacobian of the field and/or gradient.  This matrix can be inverted to determine the
    //   currents necessary to achieve a desired field and gradient.  The gradient portion of the matrix has been vectorized into 5 elements.
    //   They are [dBx/dx, dBx/dy, dBx/dz, dBy/dy, dBy/dz].
    /**
     * @brief returns the 3xN matrix mapping field at a point to the current in each of the N sources
     *  @param position is the position in the workspace the field is desired
     *
     *  This function does not check to see if the point is actually in the calibrated workspace.
     **/
    Eigen::MatrixXd fieldCurrentJacobian( const Eigen::Vector3d& position = Eigen::Vector3d::Zero() ) const;

    /**
     * @brief returns the 5xN matrix mapping the field radient at a point to the current in each of the N sources
     *  @param position is the position in the workspace the field is desired
     *
     *  The gradient matrix has been repacked, since it is symetric an has zero trace, into a five element vector.
     *    The element order is: [dBx/dx, dBx/dy, dBx/dz, dBy/dy, dBy/dz]^T
     *  This function does not check to see if the point is actually in the calibrated workspace.
     **/
    Eigen::MatrixXd gradientCurrentJacobian( const Eigen::Vector3d& position = Eigen::Vector3d::Zero() ) const;

    /**
     * @brief returns the 8xN matrix mapping of the stacked field over field gradient at a point to the current in each of the N sources
     *  @param position is the position in the workspace the field is desired
     *
     *  The gradient matrix has been repacked, since it is symetric an has zero trace, into a five element vector.
     *    The element order is: [dBx/dx, dBx/dy, dBx/dz, dBy/dy, dBy/dz]^T
     *  This function does not check to see if the point is actually in the calibrated workspace.
     **/
    Eigen::MatrixXd fieldAndGradientCurrentJacobian( const Eigen::Vector3d& position = Eigen::Vector3d::Zero() ) const;

    // This function returns the 5x3 Jacobian describing how gradient changes with position.  The gradient change is a 5x3 packing of a 3x3x3 tensor.
    //  The first column is how the gradient vector packing [dBx/dx, dBx/dy, dBx/dz, dBy/dy, dBy/dz] changes with x, the second is how it changes with y, and third is
    //  how it changes with z.
    /**
     * @brief returns the 5x3 field gradient jacobian at a point given the currents in each of the N sources
     *  @param currentVector is an orderd list of currents in each coil
     *  @param position is the position in the workspace the field is desired
     *
     *  The 3x3x3 tensor has been repacked, since the 3x3 field gradient is symetric an has zero trace, into a 5x3 element vector.
     *    The element order is:
     *     [dBx/dxdx, dBx/dxdy, dBx/dxdz
     *      dBx/dydx, dBx/dydy, dBx/dydz
     *      dBx/dzdx, dBx/dzdy, dBx/dzdz
     *      dBy/dydx, dBy/dydy, dBy/dydz
     *      dBy/dzdx, dBy/dzdy, dBy/dzdz]
     *  This function does not check to see if the point is actually in the calibrated workspace.
     **/
    Eigen::Matrix<double,5,3> gradientPositionJacobian( const Eigen::VectorXd& currentVector, const Eigen::Vector3d& position = Eigen::Vector3d::Zero() ) const;

    // This function returns the field, gradient, gradientJacobian, and field/gradient current jacobain.  It is more efficient than requesting them seporately
    /**
     * @brief calculates the field, field and field gradient jacobians, and the current jacobians at a point given the currents in each of the N sources
     *  this faunction is more efficient than calling each of the applicable functions seporately.
     *
     *  @param fieldAtPoint a pass by reference Vector to return the calculated field
     *  @param fieldGradientPositionJacobian a pass by reference matrix to return the field jacobian stacked over the field gradient jacobian
     *  @param fieldGradientCurrentJacobian a pass by reference matrix to return the field and gradient current jacobian
     *  @param currentVector is an orderd list of currents in each coil
     *  @param position is the position in the workspace the field is desired
     *
     *  The 3x3 gradient and 3x3x3 field gradient jacobian tensor has been repacked, since the 3x3 field gradient is symetric and has zero trace, into a 8x3 element vector.
     *    The element order is:
     *
     *    \f$\begin{bmatrix}
     *          \frac{\partial B_x}{\partial x} & \frac{\partial B_x}{\partial y} & \frac{\partial B_x}{\partial z} \\[0.6em]
     *          \frac{\partial B_x}{\partial y} & \frac{\partial B_y}{\partial y} & \frac{\partial B_y}{\partial z} \\[0.6em]
     *          \frac{\partial B_x}{\partial z} & \frac{\partial B_y}{\partial z} & -\left(\frac{\partial B_x}{\partial x}+\frac{\partial B_y}{\partial y}\right) \\[0.6em]
     *          \frac{\partial^2 B_x}{\partial x\partial x} & \frac{\partial^2 B_x}{\partial x\partial y} & \frac{\partial^2 B_x}{\partial x\partial z} \\[0.6em]
     *          \frac{\partial^2 B_x}{\partial y\partial x} & \frac{\partial^2 B_x}{\partial y\partial y} & \frac{\partial^2 B_x}{\partial y\partial z} \\[0.6em]
     *          \frac{\partial^2 B_x}{\partial z\partial x} & \frac{\partial^2 B_x}{\partial z\partial y} & \frac{\partial^2 B_x}{\partial z\partial z} \\[0.6em]
     *          \frac{\partial^2 B_y}{\partial y\partial x} & \frac{\partial^2 B_y}{\partial y\partial y} & \frac{\partial^2 B_y}{\partial y\partial z} \\[0.6em]
     *          \frac{\partial^2 B_y}{\partial z\partial x} & \frac{\partial^2 B_y}{\partial z\partial y} & \frac{\partial^2 B_y}{\partial z\partial z} \\[0.6em]
     *       \end{bmatrix}\f$
     *
     *  This function does not check to see if the point is actually in the calibrated workspace.
     **/
    void fullMagneticState( Eigen::Vector3d& fieldAtPoint, Eigen::Matrix<double,8,3>& fieldGradientPositionJacobian, Eigen::MatrixXd& fieldGradientCurrentJacobian, const Eigen::VectorXd& current, const Eigen::Vector3d& position = Eigen::Vector3d::Zero() ) const;

     /**
     * @brief calculates the field, field and field gradient jacobians, and the current jacobians at a point given the currents in each of the N sources
     *  this faunction is more efficient than calling each of the applicable functions seporately.
     *
     *  @param currentVector is an orderd list of currents in each coil
     *  @param position is the position in the workspace the field is desired
     *
     *  @return the information on the magnetic field at a point in the form of a MagneticState object.
     *
     */
    MagneticState fullMagneticState(const Eigen::VectorXd& currentVector, const Eigen::Vector3d& position = Eigen::Vector3d::Zero()) const;

    // These two functions convert a vector packing of the gradient to a matrix packing and vice versa
    /**
     * @brief converts a 5x1 gradient vector into a symetric zero trace 3x3 matrix
     *  @param gradVector is the desired vector to be remaped
     */
    static Eigen::Matrix3d remapGradientVector(const Vector5d& gradVector );

    /**
     * @brief converts a 3x3 gradient marix into the 5x1 gradient vector
     *  @param gradMatrix is the desired matrix to be remaped
     *
     *  This function does not check to verify the matrix is indeed symetric and zero trace
     */
    static Vector5d remapGradientMatrix(const Eigen::Matrix3d& gradMatrix);


    /**
     * @brief converts a 3x1 moment vector into a 3x5 force matrix
     *  @param moment is the magnetic moment that is packed into the force matrix
     *
     *  converts a 3x1 moment vector into a 3x5 force matrix to allow easy calculation of the magnetic force by multiplication with the gradient vector
     */
    static Eigen::MatrixXd packForceMatrix(const Eigen::Vector3d& moment);

    /**
     * @brief  loads a new calibration file
     *  @param a string pointing to the location of the yaml formated calbiration file
     */
    bool loadCalibration(std::string fileName);

    /**
     * @brief  writes a new calibration file
     *  @param a string pointing to the location for the yaml formated calbiration file
     */
    bool writeCalibration(std::string fileName) const;

    /**
     * @brief  returns the number of coils in the calibration
     */
    int getNumberOfCoils() const;

    /**
     * @brief returns the number of sources for the given coil
     */
    int getNumberOfSources( unsigned int coilNum ) const;

    /**
     * @brief returns the number of coefficients for the given source
     */
    int getNumberOfCoeffients( unsigned int coilNum, unsigned int srcNum ) const;

    /**
     * @brief returns if the calibration has a DC offset
     */
    bool hasOffset() const;

    /**
     * @brief returns the calibration name
     */
    std::string getName() const;

    // check to see if a point is within the calibrated workspace
    /**
     * @brief checks to see if a point lies within the calibrated workspace
     *  @param position the desired position to check
     */
    bool pointInWorkspace( const Eigen::Vector3d& position ) const;


    /**
     * @brief getWorkSpace returns the rectanglar extent of the calibrated workspace
     */
    MagneticWorkSpace getWorkSpace() const;

    /**
     * @brief sets the magnetic workspace to the desired specification
     */
    void setWorkSpace(const MagneticWorkSpace& ws);

    /**
     * @brief enables or disables the use of the offset, if one exists.
     */
    void useOffset(bool offsetOn );

    /**
     * @brief adds an offset to the system and enables it.
     */
    void useOffset(const ScalorPotential& newOffset );

    /**
     * @brief returns if the offset is enabled or disabled.
     */
    bool useOffset() const;

    /**
     *  \brief The calibration_constraints enum defines what constraints are active during calibration and how they are applied
     */
    enum calibration_constraints
    {
        UNIT_HEADING_ONLY, /**<  Constrains the azimuth of the potentials to be unit length */
        HEADING_AND_POSITION, /**< Constrains the azimuth of potentials to be unit length and enforces that the positions lie in a spherical anulous outside of the measured data and prevents them from going to infinity. */
        HEADING_THEN_POSITION /**< First solves with the heading constraints, then it resolves pushing any sources that lie inside the measurement data region out of a bounding circle of the data. */
    };

    /**
     * @brief calibrate  Performs a system calibration based on gathered data and a current guess as to the scalor potential structure.
     * @param calibrationName The name of the calibration.
     * @param dataList The list of magnetic measurements.
     * @param printProgress Boolean to idenify if the convergance progress should be printed with cout.
     * @param printStats Boolean to identify if the convergence statistics shoudl be printed with cout once completed.
     * @param constraint Identifies what kind of constraits should be applied to the positions and headings during convergance.
     * @param minimumSourceToCenterDistance Identifies the minimum acceptable distance between the center of the workspace and any scalor potential source.  Default is just outside workspace volume.
     * @param maximumSourceToCenterDistance Identifies the maximum acceptable distance between the center of the workspace and any scalor potential source.  Default is 10 times the initial guess distances.
     * @param converganceTolerance Specifies the convergance tolerance.  Default is 1e-12.
     * @param maxIterations  Specifies the maximum number of interations to used to prvent the loop from infinite cycles.
     * @param numberOfConvergedIterations Specifies the number of sequential converged passes to prevent a false positive in convergance.
     */
    void calibrate(std::string calibrationName, const std::vector<MagneticMeasurement>& dataList, bool printProgress = true, bool printStats = true, calibration_constraints constraint = HEADING_THEN_POSITION, double minimumSourceToCenterDistance = -1, double maximumSourceToCenterDistance = -1, double converganceTolerance = 1e-12, int maxIterations = 10000, int numberOfConvergedIterations = 1 );

    /**
     *
     *  \brief Prints to the terminal (cout) statistics describing how well the system reproduces the magnetic measurements provided.
     *  \param dataList  A set of magnetic measurements to compair the calibration to.
     *
     *  The output of this function provides the $R^2$ value, the mean error, and the square root of covariance (standard deviation) of the error in both a millitesla and normalized percentage basis.
     */
    void printStats( const std::vector< MagneticMeasurement >& dataList) const;

protected:

    ElectromagnetCalibration(); /**< Default Constructor only available to inheriting classes **/


    MagneticWorkSpace workSpace;

    std::vector<ScalorPotential> coilList;
    ScalorPotential offset;

    bool use_offset;

    std::string name;

    bool checkSourcePositions(bool printWarning = false) const;

private:

    //  Below are functions and parameters that aid in the calibration
    void applyPHI(const Eigen::VectorXd& PHI );
    void obtainPHI( Eigen::VectorXd& PHI );

    void packError( Eigen::VectorXd& error, const std::vector<MagneticMeasurement> & dataList );
    void packErrorJacobian( Eigen::MatrixXd& jacobian, const std::vector<MagneticMeasurement> & dataList );

    void linearLeastSquareCoeffFit( const std::vector<MagneticMeasurement>& dataList);

    int numberOfParameters;
    int numberOfConstraints;
    int numberOfMeasurements;
    int numberOfSources;

    int nConst;

    bool minRadIsActive;
    double rMinSq, rMaxSq, posWeight;
    Eigen::Vector3d pCenter;

};

#endif // ElectromagnetCalibration_H
