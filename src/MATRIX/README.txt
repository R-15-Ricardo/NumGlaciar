Esta es una nueva versión con menos fugas de memoria :D

Yo no revisaria "matrix.hpp". Es un archivo demasiado largo y lo unico que hace es abstraer 
todo lo feo de las matrices en c++ de una manera que se parece a <Armadillo> (Es la libreria
que a mí me gusta usar), pero como no puedo usarla, entonces hice algo "similar pero peor" ya
que parece que aún vamos a trabajar mucho con matrices.

De todas formas, esto es lo único que se necesita para saber que es lo que "matrix.hpp" hace:

---------------------------------------------------------------------------------------------

Guide for the <LINALG/matrix> library

//Enum for default types. 
enum MatrixType {
	zeros = 0, ones = 1, identity = -1, arbitrary = 2, empty = 404,
	colvector = 3, rowvector = 4
};

//Constructors
Matrix();
Matrix(unsigned int _rows, unsigned int _cols);
Matrix(unsigned int _rows, unsigned int _cols, int _type);
Matrix(const std::string& filePath);
Matrix(const std::initializer_list<std::initializer_list<double>> list);

//Memory operations
// - - - setters - - -
void setSize(unsigned int _rows, unsigned int _cols);
void setType(int _type);
bool getFromFile(const std::string& filePath);
// - - - getters - - -
unsigned int getRows() const;
unsigned int getCols() const;
int getType() const;
Matrix<type> row(unsigned int) const; //rowvec matrix
Matrix<type> col(unsigned int) const; //colvec matrix

//Operations
Matrix<type> T() const; 
void swapRows(unsigned int, unsigned int); 
void swapCols(unsigned int, unsigned int); 
void insertRow(const Matrix<type>&, unsigned int);
void insertCol(const Matrix<type>&, unsigned int);

---------------------------------------------------------------------------------------------

Todo lo demás funciona de manera relativamente intuitiva, excepto los vectores que se
acceden de maneras distintas dependiendo de si son vectores fila o vectores columna