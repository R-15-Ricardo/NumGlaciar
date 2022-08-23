#pragma once

#ifndef iostream 
	#include <iostream>
#endif

#ifndef string
	#include <string>
#endif

#define LOG(A) std::cout<<A<<std::endl

//-------------------------Class definitions-------------------------

template <class dT> // dT intended to be int float double and variants
class Matrix {
	public:
		//Matrix default types for contruction and coparisson
		enum MatrixType {
			zeros = 0, ones = 1, identity = -1, arbitrary = 2, empty = 404,
			colvector = 3, rowvector = 4
		};

		//Constructors
		Matrix() {
			// LOG("Matrix created");
		}
		Matrix(unsigned int _rows, unsigned int _cols);
		Matrix(unsigned int _rows, unsigned int _cols, int _type);
		Matrix(const std::string& filePath);
		Matrix(const std::initializer_list<std::initializer_list<double>> list);
		Matrix(const Matrix<dT>& src); //Copy constructor

		//Destructor
		~Matrix();

		//Getters
		unsigned int getRows() const {
			return this->rows;
		}

		unsigned int getCols() const {
			return this->cols;
		}

		int getType() const {
			return this->type;
		}
		Matrix<dT> row(unsigned int at) const; 
		Matrix<dT> col(unsigned int at) const;

		//Public allocators
		void setSize(unsigned int _rows, unsigned int _cols);
		void setType(int _type);
		bool getFromFile(const std::string& filePath);

		//Operators override
		Matrix<dT> operator+(const Matrix<dT>& B) const; 
		Matrix<dT> operator-(const Matrix<dT>& B) const; 
		Matrix<dT> operator*(const Matrix<dT>& B) const; 
		Matrix<dT>& operator=(const Matrix<dT>& src); //copy assigment (no idea of what im doing)
		dT* operator[](const unsigned int at) const;

		//Matrix operations
		Matrix<dT> T() const; 
		void swapRows(unsigned int i, unsigned int j); 
		void swapCols(unsigned int i, unsigned int j); 
		void insertRow(const Matrix<dT>& newRow, unsigned int at);
		void insertCol(const Matrix<dT>& newCol, unsigned int at);

		template <class U>
		friend std::ostream& operator<<(std::ostream& os, const Matrix<U>& mat);

		template <class U>
		friend Matrix<U> operator*(int, const Matrix<U>&);

		template <class U>
		friend Matrix<U> operator*(float, const Matrix<U>&);

		template <class U>
		friend Matrix<U> operator*(double, const Matrix<U>&);


	private:
		unsigned int rows {0};
		unsigned int cols {0};
		int type {MatrixType::empty};
		dT** mem {nullptr};

		//Memory private allocators
		dT** allocateMem(unsigned int _rows, unsigned int _cols) {
			dT** newMem = new dT* [_rows];
			for (int i = 0; i < _rows; i++)
				newMem[i] = new dT [_cols];
			//LOG("Memory Allocated "<<(_rows*sizeof(dT)*_cols*sizeof(dT))/1000<<"kb");
			
			return newMem;
		}

		void freeMem() { 
			for (int i = 0; i < this->rows; i++)
				delete[] this->mem[i];

			delete[] this->mem;

			//LOG("Memory freed "<<(this->rows*sizeof(dT)*this->cols*sizeof(dT))/1000<<"kb");
		}

		//Type private setters
		void fillOnes() {
			for (int i = 0; i <this->rows; i++)
				for (int j = 0; j < this->cols; j++)
					this->mem[i][j] = 1;
		}

		void fillZeros() {
			for (int i = 0; i <this->rows; i++)
				for (int j = 0; j < this->cols; j++)
					this->mem[i][j] = 0;

		}
		void fillIdentity() {
			for (int i = 0; i <this->rows; i++)
				for (int j = 0; j < this->cols; j++)
					this->mem[i][j] = 0;

			int n = (this->rows < this->cols) ? this->rows : this->cols;

			for (int i = 0; i < n; i++)
				this->mem[i][i] = 1;
		}
};

typedef Matrix<double> dMatrix;
typedef Matrix<double> dVector;

//-------------------------Class definitions-------------------------

// - - - - - - - - Constructors - - - - - - - -
template <class dT>
Matrix<dT>::Matrix(unsigned int _rows, unsigned int _cols) : rows(_rows), cols(_cols)
{
	// LOG("Matrix created");

	this->type = MatrixType::arbitrary;

	this->mem = allocateMem(_rows, _cols);
	if (this->mem == nullptr)
		perror("Matrix(): Out of memory!");
}

template <class dT>
Matrix<dT>::Matrix(unsigned int _rows, unsigned int _cols, int _type) : rows(_rows), cols(_cols)
{
	// LOG("Matrix created");

	this->mem = allocateMem(_rows, _cols);
	if (this->mem == nullptr)
		perror("Matrix(): Out of memory!");

	switch (_type)
	{
		case MatrixType::ones:
			this->fillOnes();
			break;

		case MatrixType::zeros:
			this->fillZeros();
			break;

		case MatrixType::identity:
			this->fillIdentity();
			break;
		
		default: 
			break;
	}

	this->type = _type;
}

template <class dT>
Matrix<dT>::Matrix(const std::string& filePath)
{
	// LOG("Matrix created");
	this->getFromFile(filePath);
}

template <class dT>
Matrix<dT>::Matrix(const std::initializer_list<std::initializer_list<double>> L)
{
	// LOG("Matrix created");

	//Only this constructor is independent of Matrix<type>allocateMem()
	this->rows = L.size();
	this->mem = new dT* [this->rows];

	int i = 0;
	int prevM = 0, nowM = 0;
	for (auto l: L)
	{
		if (nowM != prevM)
			perror("Matrix(): List size does not match");

		int j = 0;
		nowM = l.size();
		this->mem[i] = new dT [nowM];
		for (auto x: l)
			this->mem[i][j++] = static_cast<dT>(x);
		i++;
		prevM = nowM;
	}
	this->cols = nowM;
	this->type = MatrixType::arbitrary;
}

template <class dT>
Matrix<dT>::Matrix(const Matrix<dT>& src)
{
	// LOG("Matrix created");

	this->rows = src.getRows();
	this->cols = src.getCols();
	this->type = src.getType();

	this->mem = this->allocateMem(this->rows, this->cols);

	for (int i = 0; i < this->rows; i++)
		for (int j = 0; j < this->cols; j++)
			this->mem[i][j] = src[i][j];
}

// - - - - - - - - Destructor - - - - - - - -

template <class dT>
Matrix<dT>::~Matrix()
{
	// LOG("Matrix deleted");
	this->freeMem();
}

// - - - - - - - - Getters - - - - - - - -

template <class dT>
Matrix<dT> Matrix<dT>::row(unsigned int at) const
{
	Matrix<dT> out(1,this->cols,MatrixType::rowvector);
	for (int i = 0; i < this->cols; i++)
		out[0][i] = this->mem[at][i]; 

	return out;
}

template <class dT>
Matrix<dT> Matrix<dT>::col(unsigned int at) const
{
	Matrix<dT> out(this->rows,1,MatrixType::colvector);
	for (int i = 0; i < this->rows; i++)
		out[i][0] = this->mem[i][at]; 
	
	return out;
}

// - - - - - - - - Public allocators - - - - - - - -
template <class dT>
void Matrix<dT>::setSize(unsigned int _rows, unsigned int _cols)
{
	return;
}

template <class dT>
void Matrix<dT>::setType(int _type)
{
	if (this->type != MatrixType::empty)
	{
		switch (_type)
		{
			//Basic setters (size independent)
			case MatrixType::identity:
				this->fillIdentity();
				break;

			case MatrixType::ones:
				this->fillOnes();
				break;
			
			case MatrixType::zeros:
				this->fillZeros();
				break;

			//Morphing setters
			case MatrixType::rowvector:
				if (this->rows > 1)
				{
					this->freeMem();
					this->mem = this->allocateMem(1, this->cols);
					this->rows = 1; 
				}
				break;

			case MatrixType::colvector:
				if (this->cols > 1)
				{
					this->freeMem();
					this->mem = this->allocateMem(this->rows, 1);
					this->cols = 1; 
				}
				break;

			case MatrixType::empty:
				this->freeMem();
				this->rows = this->cols = 0;

			default:
				perror("setType(): Matrix type not defined!");
				break;
		}

		this->type = _type;
	}
	else if (this->type == MatrixType::empty && _type == MatrixType::empty)
		return;
	else
		perror("setType(): Can't assing values to empty matrix!");

}

template <class dT>
bool Matrix<dT>::getFromFile(const std::string& filePath)
{
	this->rows = this->cols = 0;
	this->type = MatrixType::empty;

	if (this->mem != nullptr)
		this->freeMem();

	FILE *f = fopen(filePath.c_str(), "rb");

	if (!f) 
	{
		perror("getFromFile(): Failed to load file");
		return false;
	}

	//This makes the method slightly slower, but much prettier
	std::string inFileType = filePath.substr(filePath.size()-3,filePath.size());

	if (inFileType == "mat")
		this->type = MatrixType::arbitrary;
	else if (inFileType == "vec")
		this->type = MatrixType::colvector;
	else
	{
		perror("getFromFile(): File extension not supported!");
		return false;
	}

	fread(&this->rows, sizeof(int), 1, f);

	if (this->type == MatrixType::arbitrary)
    	fread(&this->cols, sizeof(int), 1, f);
	else 
		this->cols = 1;

	//Straight of the last version, so still completely unintelligible
	this->mem = new dT* [this->rows];
	for (int i = 0; i<this->rows; i++)
	{
		this->mem[i] = new dT [this->cols];
		fread(this->mem[i], sizeof(dT), (this->cols), f);
	}

	fclose(f);
	return true;
}


// - - - - - - - - Operators overload - - - - - - - -

template <class dT>
Matrix<dT> Matrix<dT>::operator+(const Matrix<dT>& B) const
{
	if (this->rows != B.getRows() || this->cols != B.getCols())
	{
		perror("operator+(): Matrix size doesn't match!");
		return Matrix();
	}

	int _type;
	
	if (this->type == MatrixType::rowvector || this->type == MatrixType::colvector)
		_type = this->type;
	else if (this->type == MatrixType::zeros || this->type == MatrixType::zeros)
		_type = this->type + B.getType(); 
	else 
		_type = MatrixType::arbitrary;

	Matrix<dT> retMat(this->rows, this->cols, _type);
	for (int i = 0; i < this->rows; i++)
		for (int j = 0; j < this->cols; j++)
			retMat[i][j] = this->mem[i][j] + B[i][j];
	
	return retMat;
}

template <class dT>
Matrix<dT> Matrix<dT>::operator-(const Matrix<dT>& B) const
{
	if (this->rows != B.getRows() || this->cols != B.getCols())
	{
		perror("operator-(): Matrix size doesn't match!");
		return Matrix();
	}

	int _type;
	if (this->type != MatrixType::arbitrary && 
		this->type != MatrixType::colvector &&
		this->type != MatrixType::rowvector && 
		this->type == B.getType())
		_type = MatrixType::zeros;
	
	else 
		_type = this->type;

	Matrix<dT> retMat(this->rows, this->cols, _type);
	for (int i = 0; i < this->rows; i++)
		for (int j = 0; j < this->cols; j++)
			retMat[i][j] = this->mem[i][j] - B[i][j];

	return retMat;
}

template <class dT>
Matrix<dT> Matrix<dT>::operator*(const Matrix<dT>& B) const
{
	if (this->cols != B.getRows())
	{
		perror("operator*(): Matrix size doesn't match");
		return Matrix<dT>();
	}
	
	int _type;
	if (B.getType() == MatrixType::colvector)
		_type = MatrixType::colvector;

	else _type = this->type;

	Matrix<dT> retMat(this->rows, B.getCols(), _type);

	for (int i = 0; i < this->rows; i++)
	{
		for (int j = 0; j < B.getCols(); j++)
		{
			dT sum=0;
			for (int k = 0; k < this->cols; k++)
				sum += this->mem[i][k]*B[k][j];

			retMat[i][j] = sum;
		}
	}

	return retMat;
}

template <class dT>
Matrix<dT>& Matrix<dT>::operator=(const Matrix<dT>& src) 
{
	// LOG("operator=(): Copy assigment called!");

	if (this->mem != nullptr)
		this->freeMem();
	
	this->rows = src.getRows();
	this->cols = src.getCols();
	this->type = src.getType();
	this->mem = this->allocateMem(src.getRows(), src.getCols());

	for (int i = 0; i < this->rows; i++)
		for (int j = 0; j < this->cols; j++)
			this->mem[i][j] = src[i][j];

	return *this;
}

template <class dT>
dT* Matrix<dT>::operator[](const unsigned int at) const 
{
	return this->mem[at];
}

// - - - - - - - - Matrix operations - - - - - - - -

template <class dT>
Matrix<dT> Matrix<dT>::T() const
{
	int _type;
	switch (this->type)
	{
	case MatrixType::ones:
		_type = this->type;
		break;
	case MatrixType::zeros:
		_type = this->type;
		break;
	case MatrixType::rowvector:
		_type = MatrixType::colvector;
		break;
	case MatrixType::colvector:
		_type = MatrixType::rowvector;
		break;
	default:
		_type = MatrixType::arbitrary;
		break;
	}

	Matrix<dT> transposed(this->cols, this->rows, _type);
	for (int i = 0; i < this->rows; i++)
	{
		for (int j = 0; j < this->cols; j++)
			transposed[j][i] = this->mem[i][j];
	}

	return transposed;
}



template <class dT>
void Matrix<dT>::swapRows(unsigned int i, unsigned int j)
{
	dT aux;
	for (int k = 0; k < this->cols; k++)
	{
		aux = this->mem[i][k];
		this->mem[i][k] = this->mem[j][k];
		this->mem[j][k] = aux;
	}
}

template <class dT>
void Matrix<dT>::swapCols(unsigned int i, unsigned int j)
{
	dT aux;
	for (int k = 0; k < this->cols; k++)
	{
		aux = this->mem[k][i];
		this->mem[k][i] = this->mem[k][j];
		this->mem[k][j] = aux;
	}
}

template <class dT>
void Matrix<dT>::insertRow(const Matrix<dT>& newRow, unsigned int at)
{
	if (newRow.getType() != MatrixType::rowvector || this->cols != newRow.getCols())
	{
		perror("insertRow(): Vector type doesn't match");
		return;
	}

	for (int i=0; i < this->cols; i++)
		this->mem[at][i] = newRow[0][i];
}

template <class dT>
void Matrix<dT>::insertCol(const Matrix<dT>& newCol, unsigned int at)
{
	if (newCol.getType() != MatrixType::colvector || this->rows != newCol.getRows())
	{
		perror("insertCol(): Vector type doesn't match");
		return;
	}

	for (int i=0; i < this->rows; i++)
		this->mem[i][at] = newCol[i][0];
}

// - - - - - - - - operator overload (friends) - - - - - - - -

template <class dT>
std::ostream& operator<<(std::ostream& stream, const Matrix<dT>& a)
{
	if (a.type == Matrix<dT>::MatrixType::empty)
	{	
		stream<<"[/]";
		return stream;
	}
	//this is terrible. oh well
	if (a.rows <= 8 && a.cols <= 8)
	{
		for (int i = 0; i < a.rows; i++)
		{
			stream<<"| ";
			for (int j = 0; j < a.cols; j++) 
			{
				stream<<a.mem[i][j];
				stream<<" ";
			}
			stream<<"|\n";
		}
	}

	else if (a.rows <= 8 && a.cols > 8)
	{
		for (int i = 0; i < a.rows; i++)
		{
			stream<<"| ";
			for (int j = 0; j < 4; j++) 
			{
				stream<<a.mem[i][j];
				stream<<" ";
			}
			stream<<"... ";
			for (int j = a.cols-4; j < a.cols; j++) 
			{
				stream<<a.mem[i][j];
				stream<<" ";  
			} 
			stream<<"|\n";
		}
	}

	else if (a.rows > 8 && a.cols <= 8)
	{
		for (int i = 0; i < 4; i++)
		{
			stream<<"| ";
			for (int j = 0; j < a.cols; j++) 
			{
				stream<<a.mem[i][j];
				stream<<" ";
			}
			stream<<"|\n";
		}
		stream<<"    .\n    .\n    .\n";
		for (int i = a.rows-4; i < a.rows; i++)
		{
			stream<<"| ";
			for (int j = 0; j < a.cols; j++) 
			{
				stream<<a.mem[i][j];
				stream<<" ";
			}
			stream<<"|\n";
		}
	}

	else
	{
		for (int i = 0; i < 4; i++)
		{
			stream<<"| ";
			for (int j = 0; j < 4; j++) 
			{
				stream<<a.mem[i][j];
				stream<<" "; 
			}
			stream<<"... ";
			for (int j = a.cols-4; j < a.cols; j++) 
			{
				stream<<a.mem[i][j];
				stream<<" ";
			} 
			stream<<"|\n";
		}
		stream<<"    .\n    .\n    .\n";
		for (int i = a.rows-4; i < a.rows; i++)
		{
			stream<<"| ";
			for (int j = 0; j < 4; j++) 
			{
				stream<<a.mem[i][j];
				stream<<" ";   
			}
			stream<<"... ";
			for (int j = a.cols-4; j < a.cols; j++) 
			{
				stream<<a.mem[i][j];
				stream<<" "; 
			} 
			stream<<"|\n";
		}
	}
	

	return stream;
}

template <class dT>
Matrix<dT> operator*(int escalar, const Matrix<dT>& a)
{
	int _type;
	if (escalar == 0)
		_type = Matrix<dT>::MatrixType::zeros;
	else
		_type = a.type;

	Matrix<dT> retMat(a.rows, a.cols,_type);
	for (int i = 0; i < a.rows; i++)
	{
		for (int j = 0; j < a.cols; j++)
		{
			retMat[i][j] = escalar * a.mem[i][j];
		}
	}

	return retMat;
}

template <class dT>
Matrix<dT> operator*(float escalar, const Matrix<dT>& a)
{
	int _type;
	if (escalar == 0)
		_type = Matrix<dT>::MatrixType::zeros;
	else
		_type = a.type;

	Matrix<dT> retMat(a.rows, a.cols,_type);
	for (int i = 0; i < a.rows; i++)
	{
		for (int j = 0; j < a.cols; j++)
		{
			retMat[i][j] = escalar * a.mem[i][j];
		}
	}

	return retMat;
}

template <class dT>
Matrix<dT> operator*(double escalar, const Matrix<dT>& a)
{
	int _type;
	if (escalar == 0)
		_type = Matrix<dT>::MatrixType::zeros;
	else
		_type = a.type;

	Matrix<dT> retMat(a.rows, a.cols,_type);
	for (int i = 0; i < a.rows; i++)
	{
		for (int j = 0; j < a.cols; j++)
		{
			retMat[i][j] = escalar * a.mem[i][j];
		}
	}

	return retMat;
}
