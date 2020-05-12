#include <mpi.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <thread>
#include <chrono>

using namespace std;

struct Offset {
    int I, J, K;
};

int calculateFlatIndex(int pointsPerEdge, int i, int j, int k) {
    return i * pointsPerEdge * pointsPerEdge + j * pointsPerEdge + k;
}

int calculateFlatOffset(int pointsPerEdge, Offset myOffset) {
    return calculateFlatIndex(pointsPerEdge, myOffset.I, myOffset.J, myOffset.K);
}

void fillWithZeros(double* cube, int pointsPerEdge) {
    for (int i = 0; i < pointsPerEdge; i++) {
        for (int j = 0; j < pointsPerEdge; j++) {
            for (int k = 0; k < pointsPerEdge; k++) {
                cube[calculateFlatIndex(pointsPerEdge, i, j, k)] = 0;
            }
        }
    }
}

void fillByEnumerating(double* cube, int pointsPerEdge, int offset) {
    for (int i = 0; i < pointsPerEdge; i++) {
        for (int j = 0; j < pointsPerEdge; j++) {
            for (int k = 0; k < pointsPerEdge; k++) {
                int flatIndex = calculateFlatIndex(pointsPerEdge, i, j, k);
                cube[flatIndex] = offset + flatIndex;
            }
        }
    }
}

void print(double* cube, int pointsPerEdge) {
    for (int i = 0; i < pointsPerEdge; i++) {
        for (int j = 0; j < pointsPerEdge; j++) {
            for (int k = 0; k < pointsPerEdge; k++) {
                cout << cube[calculateFlatIndex(pointsPerEdge, i, j, k)] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}

bool isWholeNumber(double num)
{
    return num == static_cast<int>(num);
}

void calculateMyOffset(int myRank, int cubeDimension, int subcubeDimension, Offset &myOffset) {

    myOffset.I = 0;
    myOffset.J = 0;
    myOffset.K = 0;

    for (int i = 0; i < myRank; i++) {

        myOffset.K += subcubeDimension;

        if (myOffset.K == cubeDimension) {

            myOffset.K = 0;

            myOffset.J += subcubeDimension;

            if (myOffset.J == cubeDimension) {

                myOffset.J = 0;

                myOffset.I += subcubeDimension;
            }

        }

    }

}

void printBuffer(double* buffer, int bufferSize) {
    for (int i = 0; i < bufferSize; i++) {
        cout << buffer[i] << " ";
    }
    cout << endl;
}

void fillBufferWithZeros(double* buffer, int bufferSize) {
    for (int i = 0; i < bufferSize; i++) {
        buffer[i] = 0;
    }
}

void calculateCurrentOffset(Offset& myOffset, int i, int j, int k, Offset& currentOffset) {
    currentOffset.I = myOffset.I + i;
    currentOffset.J = myOffset.J + j;
    currentOffset.K = myOffset.K + k;

}

bool isBoundryPoint(Offset& offset, int pointsPerEdge) {

    int minIndex = 0;
    int maxIndex = pointsPerEdge - 1;
    
    if (offset.I == minIndex || offset.I == maxIndex) {
        return true;
    }

    if (offset.J == minIndex || offset.J == maxIndex) {
        return true;
    }

    if (offset.K == minIndex || offset.K == maxIndex) {
        return true;
    }

    return false;
}

double u(double x, double y, double z) {
    return x * y * z;
}

double un(int n, int i, int j, int k) {
    double x = (1.0 * i) / n;
    double y = (1.0 * j) / n;
    double z = (1.0 * k) / n;
    return u(x, y, z);
}

double f(double x, double y, double z) {
    return 1.0 + 1.0 + 1.0;
}

double fn(int n, int i, int j, int k) {
    double x = (1.0 * i) / n;
    double y = (1.0 * j) / n;
    double z = (1.0 * k )/ n;
    return f(x, y, z);
}

int main(int argc, char* argv[]) {

    // MPI-Initialization -----------------------------------------------------------------------------

    int numberOfTasks, myRank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfTasks);

    // Points per Edge Calculations -------------------------------------------------------------------------------------

    double subcubePerEdge = cbrt(numberOfTasks);

    if (!isWholeNumber(subcubePerEdge)) {
        cout << "cube root of number of tasks is not a whole number!" << endl;
        return 0;
    }

    int n = stoi(argv[1]);

    int pointsPerEdge = n + 1;

    double _pointsPerSubEdge = pointsPerEdge / subcubePerEdge;

    if (!isWholeNumber(_pointsPerSubEdge)) {
        cout << "points per subcube is not a whole number!" << endl;
        return 0;
    }

    int pointsPerSubEdge = (int) _pointsPerSubEdge;

    int pointsPerGhostEdge = 1 + pointsPerSubEdge + 1;

    // Initialize Subcube -------------------------------------------------------------

    double *subcube = new double[pointsPerSubEdge * pointsPerSubEdge * pointsPerSubEdge];

    fillWithZeros(subcube, pointsPerSubEdge);

    double* ghostCube = new double[pointsPerGhostEdge * pointsPerGhostEdge * pointsPerGhostEdge];

    fillWithZeros(ghostCube, pointsPerGhostEdge);

    // Calculate Offsets -------------------------------------------------------------------------------------

    struct Offset myOffset;

    calculateMyOffset(myRank, pointsPerEdge, pointsPerSubEdge, myOffset);

    int offsetFlatIndex = calculateFlatOffset(pointsPerEdge, myOffset);

    Offset myNormalizedOffset = {myOffset.I / pointsPerSubEdge, myOffset.J / pointsPerSubEdge, myOffset.K / pointsPerSubEdge };

    // Create Neighbours Cube -----------------------------------------------------------------------------------

    int*** neighbours  = new int** [subcubePerEdge];

    for (int i = 0; i < subcubePerEdge; i++) {
        
        neighbours[i] = new int* [subcubePerEdge];

        for (int j = 0; j < subcubePerEdge; j++) {
           
            neighbours[i][j] = new int[subcubePerEdge];
       
        }

    }
      
    for (int i = 0; i < subcubePerEdge; i++) {
        for (int j = 0; j < subcubePerEdge; j++) {
            for (int k = 0; k < subcubePerEdge; k++) {
                neighbours[i][j][k] = calculateFlatIndex(subcubePerEdge, i, j, k);
            }
        }
    }

    // Find Neighbours ---------------------------------------------------------------------------------------

    int ME = myRank;

    int UP_NEIGHBOUR = -1;
    int DOWN_NEIGHBOUR = -1;
    
    int LEFT_NEIGHBOUR = -1;
    int RIGHT_NEIGHBOUR = -1;

    int FRONT_NEIGHBOUR = - 1;
    int BACK_NEIGHBOUR = -1;

    if (myNormalizedOffset.I - 1 >= 0) {
        UP_NEIGHBOUR = neighbours[myNormalizedOffset.I - 1][myNormalizedOffset.J][myNormalizedOffset.K];
    }

    if (myNormalizedOffset.I + 1 < subcubePerEdge) {
        DOWN_NEIGHBOUR = neighbours[myNormalizedOffset.I + 1][myNormalizedOffset.J][myNormalizedOffset.K];
    }

    if (myNormalizedOffset.J - 1 >= 0) {
        LEFT_NEIGHBOUR = neighbours[myNormalizedOffset.I][myNormalizedOffset.J - 1][myNormalizedOffset.K];
    }

    if (myNormalizedOffset.J + 1 < subcubePerEdge) {
        RIGHT_NEIGHBOUR = neighbours[myNormalizedOffset.I][myNormalizedOffset.J + 1][myNormalizedOffset.K];
    }

    if (myNormalizedOffset.K - 1 >= 0) {
        BACK_NEIGHBOUR = neighbours[myNormalizedOffset.I][myNormalizedOffset.J][myNormalizedOffset.K - 1];
    }

    if (myNormalizedOffset.K + 1 < subcubePerEdge) {
        FRONT_NEIGHBOUR = neighbours[myNormalizedOffset.I][myNormalizedOffset.J][myNormalizedOffset.K + 1];
    }

    // Buffers --------------------------------------------------------------------------------------------

    int bufferSize = pointsPerSubEdge * pointsPerSubEdge;

    double* upBuffer = new double[bufferSize];
    fillBufferWithZeros(upBuffer, bufferSize);

    double* downBuffer = new double[bufferSize];
    fillBufferWithZeros(downBuffer, bufferSize);

    double* frontBuffer = new double[bufferSize];
    fillBufferWithZeros(frontBuffer, bufferSize);

    double* backBuffer = new double[bufferSize];
    fillBufferWithZeros(backBuffer, bufferSize);

    double* leftBuffer = new double[bufferSize];
    fillBufferWithZeros(leftBuffer, bufferSize);

    double* rightBuffer = new double[bufferSize];
    fillBufferWithZeros(rightBuffer, bufferSize);

    // Layer Calculations ----------------------------------------------------------------------------------

    MPI_Datatype XZLayerDataType; // For Up & Down
    MPI_Type_vector(1, pointsPerSubEdge * pointsPerSubEdge, 0, MPI_DOUBLE, &XZLayerDataType);
    MPI_Type_commit(&XZLayerDataType); // MPI_Type_vector (count,blocklength,stride,oldtype,&newtype)

    MPI_Datatype XYLayerDataType; // For Front & Back
    MPI_Type_vector(pointsPerSubEdge * pointsPerSubEdge, 1, pointsPerSubEdge, MPI_DOUBLE, &XYLayerDataType);
    MPI_Type_commit(&XYLayerDataType); // MPI_Type_vector (count,blocklength,stride,oldtype,&newtype)

    MPI_Datatype YZLayerDataType; // For Left & Right
    MPI_Type_vector(pointsPerSubEdge, pointsPerSubEdge, pointsPerSubEdge * pointsPerSubEdge, MPI_DOUBLE, &YZLayerDataType);
    MPI_Type_commit(&YZLayerDataType); // MPI_Type_vector (count,blocklength,stride,oldtype,&newtype)

    // Tag Constants --------------------------------------------------------------------------------------

    const int UP_LAYER_TAG = 1;
    const int DOWN_LAYER_TAG = 2;
    const int FRONT_LAYER_TAG = 3;
    const int BACK_LAYER_TAG = 4;
    const int LEFT_LAYER_TAG = 5;
    const int RIGHT_LAYER_TAG = 6;

    // Calculate Boundary Points

    Offset currentOffset;

    for (int i = 0; i < pointsPerSubEdge; i++) {
        for (int j = 0; j < pointsPerSubEdge; j++) {
            for (int k = 0; k < pointsPerSubEdge; k++) {
                calculateCurrentOffset(myOffset, i, j, k, currentOffset);
                if (!isBoundryPoint(currentOffset, pointsPerEdge)) {
                     subcube[calculateFlatIndex(pointsPerSubEdge, i, j, k)] = un(n, currentOffset.I, currentOffset.J, currentOffset.K);
                }
            }
        }
    }

    // Send --------------------------------------------------------------------------------

    // Test Purposes Only: fillByEnumerating(subcube, pointsPerSubEdge, myRank * pointsPerSubEdge * pointsPerSubEdge * pointsPerSubEdge);
    
    // Send Up Layer
    if (UP_NEIGHBOUR != -1) {
        MPI_Send(&subcube[0], 1, XZLayerDataType, UP_NEIGHBOUR, UP_LAYER_TAG, MPI_COMM_WORLD);
    }

    // Send Down Layer
    if (DOWN_NEIGHBOUR != -1) {
        MPI_Send(&subcube[pointsPerSubEdge * pointsPerSubEdge * (pointsPerSubEdge - 1)], 1, XZLayerDataType, DOWN_NEIGHBOUR, DOWN_LAYER_TAG, MPI_COMM_WORLD);
    }

    // Send Front Layer
    if (FRONT_NEIGHBOUR != -1) {
        MPI_Send(&subcube[pointsPerSubEdge - 1], 1, XYLayerDataType, FRONT_NEIGHBOUR, FRONT_LAYER_TAG, MPI_COMM_WORLD);
    }

    // Send Back Layer
    if (BACK_NEIGHBOUR != -1) {
        MPI_Send(&subcube[0], 1, XYLayerDataType, BACK_NEIGHBOUR, BACK_LAYER_TAG, MPI_COMM_WORLD);
    }

    // Send Left Layer
    if (LEFT_NEIGHBOUR != -1) {
        MPI_Send(&subcube[0], 1, YZLayerDataType, LEFT_NEIGHBOUR, LEFT_LAYER_TAG, MPI_COMM_WORLD);
    }

    // Send Right Layer
    if (RIGHT_NEIGHBOUR != -1) {
        MPI_Send(&subcube[pointsPerSubEdge * pointsPerSubEdge - pointsPerSubEdge], 1, YZLayerDataType, RIGHT_NEIGHBOUR, RIGHT_LAYER_TAG, MPI_COMM_WORLD);
    }
 
    // Receive --------------------------------------------------------------------------------

    MPI_Status downBufferStatus;

    // Receive Down Layer
    if (DOWN_NEIGHBOUR != -1) {
        MPI_Recv(downBuffer, pointsPerSubEdge * pointsPerSubEdge, MPI_DOUBLE, DOWN_NEIGHBOUR, UP_LAYER_TAG, MPI_COMM_WORLD, &downBufferStatus);
    }

    MPI_Status upBufferStatus;

    // Receive Up Layer
    if (UP_NEIGHBOUR != -1) {
        MPI_Recv(upBuffer, pointsPerSubEdge * pointsPerSubEdge, MPI_DOUBLE, UP_NEIGHBOUR, DOWN_LAYER_TAG, MPI_COMM_WORLD, &upBufferStatus);
    }

    MPI_Status frontBufferStatus;

    // Receive Front Layer
    if (FRONT_NEIGHBOUR != -1) {
        MPI_Recv(frontBuffer, pointsPerSubEdge * pointsPerSubEdge, MPI_DOUBLE, FRONT_NEIGHBOUR, BACK_LAYER_TAG, MPI_COMM_WORLD, &frontBufferStatus);
    }

    MPI_Status backBufferStatus;

    // Receive Back Layer
    if (BACK_NEIGHBOUR != -1) {
        MPI_Recv(backBuffer, pointsPerSubEdge * pointsPerSubEdge, MPI_DOUBLE, BACK_NEIGHBOUR, FRONT_LAYER_TAG, MPI_COMM_WORLD, &backBufferStatus);
    }

    MPI_Status leftBufferStatus;

    // Receive Left Layer
    if (LEFT_NEIGHBOUR != -1) {
        MPI_Recv(leftBuffer, pointsPerSubEdge * pointsPerSubEdge, MPI_DOUBLE, LEFT_NEIGHBOUR, RIGHT_LAYER_TAG, MPI_COMM_WORLD, &leftBufferStatus);
    }

    MPI_Status rightBufferStatus;

    // Receive Left Layer
    if (RIGHT_NEIGHBOUR != -1) {
        MPI_Recv(rightBuffer, pointsPerSubEdge * pointsPerSubEdge, MPI_DOUBLE, RIGHT_NEIGHBOUR, LEFT_LAYER_TAG, MPI_COMM_WORLD, &rightBufferStatus);
    }
    
    // DEBUG INFO -------------------------------------------------------------------------------------------

    std::this_thread::sleep_for(std::chrono::milliseconds(myRank * 250));

    if (myRank == 0) {
    
        cout << "-------- DEBUG INFO --------" << endl;

        cout << endl;
        cout << "Number of Tasks: " << numberOfTasks << endl;
        cout << "Subcube Per Edge: " << subcubePerEdge << endl;
        cout << "n: " << n << endl;
        cout << "Points per Edge: " << pointsPerEdge << endl;
        cout << "Points per Sub Edge: " << pointsPerSubEdge << endl;
        cout << "Points per Ghost Edge: " << pointsPerGhostEdge << endl;
        cout << endl;

        // Test Purposes
  
    
    }

    cout << "My Rank: " << myRank << endl;

    cout << "My Offset: " << myOffset.I << " " << myOffset.J << " " << myOffset.K << endl;

    cout << "My Flat Offset : " << offsetFlatIndex << endl;

    cout << "My Normalized Offset: " << myNormalizedOffset.I << " " << myNormalizedOffset.J << " " << myNormalizedOffset.K << endl;

    cout << "Neighbours:" << " U " << UP_NEIGHBOUR << " D " << DOWN_NEIGHBOUR << " L " << LEFT_NEIGHBOUR << " R " << RIGHT_NEIGHBOUR << " F " << FRONT_NEIGHBOUR << " B " << BACK_NEIGHBOUR << endl;

    cout << "Up Buffer: " << endl;

    printBuffer(upBuffer, bufferSize);

    cout << "Down Buffer: " << endl;

    printBuffer(downBuffer, bufferSize);

    cout << "Front Buffer: " << endl;

    printBuffer(frontBuffer, bufferSize);

    cout << "Back Buffer: " << endl;

    printBuffer(backBuffer, bufferSize);

    cout << "Left Buffer: " << endl;

    printBuffer(leftBuffer, bufferSize);

    cout << "Right Buffer: " << endl;

    printBuffer(rightBuffer, bufferSize);

    cout << endl;

    print(subcube, pointsPerSubEdge);

    //print(ghostCube, pointsPerGhostEdge);

    // Finalize -------------------------------------------------------------------------------------------

    MPI_Type_free(&XYLayerDataType);
    MPI_Type_free(&XZLayerDataType);
    MPI_Type_free(&YZLayerDataType);

    MPI_Finalize();

    return 0;
}