// supported Type:     float, double, complex<float>, complex<double>      

#include "allocate.h"
#include <iostream>
#include <math.h>
#include <cstring>
using namespace std;
 
template <class Type>
class matrix;
template <class Type>
class vector {
    int n = 0;
    Type* x = NULL; 
  public:
    vector() {
        //printf("vector construct()\n");
        n=0;
        if(x != NULL) { free(x); x = NULL; }
    }
    vector(int m) {
        int i,j;
#if AVX512F > 0
        typedef typename std::conditional<std::is_same_v<double,Type> || std::is_same_v<complex<double>,Type>,__m512d,__m512>::type avxtype;
        alignas(ALIGN) avxtype mw01_vecmat_a;
#elif AVX > 0
        typedef typename std::conditional<std::is_same_v<double,Type> || std::is_same_v<complex<double>,Type>,__m256d,__m256>::type avxtype;
        alignas(ALIGN) avxtype mw01_vecmat_a;
#endif
        //printf("vector construct(m)\n");
        n=m;
        if(x != NULL) { free(x); x = NULL; }

        x = allocate1D<Type>(n);
#if !defined(AVX) || AVX == 0
        memset(x,0,n*sizeof(Type));
#elif AVX512F == 0
        j = 32/sizeof(Type);
        if constexpr(std::is_same_v<double,Type> || std::is_same_v<complex<double>,Type>) {
            mw01_vecmat_a = _mm256_setzero_pd();
            for(i=0;i<n;i+=j) _mm256_store_pd((double *)&x[i],mw01_vecmat_a);
        } else if constexpr(std::is_same_v<float,Type> || std::is_same_v<complex<float>,Type>) {
            mw01_vecmat_a = _mm256_setzero_ps();
            for(i=0;i<n;i+=j) _mm256_store_ps((float *)&x[i],mw01_vecmat_a);
        }
#else
        j = 64/sizeof(Type);
        if constexpr(std::is_same_v<double,Type> || std::is_same_v<complex<double>,Type>) {
            mw01_vecmat_a = _mm512_setzero_pd();
            for(i=0;i<n;i+=j) _mm512_store_pd((double *)&x[i],mw01_vecmat_a);
        } else if constexpr(std::is_same_v<float,Type> || std::is_same_v<complex<float>,Type>) {
            mw01_vecmat_a = _mm512_setzero_ps();
            for(i=0;i<n;i+=j) _mm512_store_ps((float *)&x[i],mw01_vecmat_a);
        }
#endif
        //for(int i=0;i<n;i++) x[i] = Type(0.);
    }
    vector(int m,Type* y)                                 // y must also be aligned array
    {
        int i,j;
        //printf("vector construct(m,y)\n");
        n=m;
        if(x != NULL) { free(x); x = NULL; }
        x = allocate1D<Type>(n);
#if !defined(AVX) || AVX == 0
        memcpy(x,y,n*sizeof(Type));
#elif AVX512F == 0
        j = 32/sizeof(Type);
        if constexpr(std::is_same_v<double,Type> || std::is_same_v<complex<double>,Type>) {
            for(i=0;i<n;i+=j) _mm256_store_pd((double *)&x[i],_mm256_load_pd((double *)&y[i]));
        } else if constexpr(std::is_same_v<float,Type> || std::is_same_v<complex<float>,Type>) {
            for(i=0;i<n;i+=j) _mm256_store_ps((float *)&x[i],_mm256_load_ps((float *)&y[i]));
        }
#else
        j = 64/sizeof(Type);
        if constexpr(std::is_same_v<double,Type> || std::is_same_v<complex<double>,Type>) {
            for(i=0;i<n;i+=j) _mm512_store_pd((double *)&x[i],_mm512_load_pd((double *)&y[i]));
        } else if constexpr(std::is_same_v<float,Type> || std::is_same_v<complex<float>,Type>) {
            for(i=0;i<n;i+=j) _mm512_store_ps((float *)&x[i],_mm512_load_ps((float *)&y[i]));
        }
#endif
    }

    vector(const vector<Type>& v) {                             // copy constructor
        int i,j;                                                // class a = b;
        //printf("vector construct &\n");
        if(this == &v) return;
        n = v.n;
        if(x != NULL) { free(x); x = NULL; }
        x = allocate1D<Type>(n);
#if !defined(AVX) || AVX == 0
        memcpy(x,v.x,n*sizeof(Type));
#elif AVX512F == 0
        j = 32/sizeof(Type);
        if constexpr(std::is_same_v<double,Type> || std::is_same_v<complex<double>,Type>) {
            for(i=0;i<n;i+=j) _mm256_store_pd((double *)&x[i],_mm256_load_pd((double *)&v.x[i]));
        } else if constexpr(std::is_same_v<float,Type> || std::is_same_v<complex<float>,Type>) {
            for(i=0;i<n;i+=j) _mm256_store_ps((float *)&x[i],_mm256_load_ps((float *)&v.x[i]));
        }
#else
        j = 64/sizeof(Type);
        if constexpr(std::is_same_v<double,Type> || std::is_same_v<complex<double>,Type>) {
            for(i=0;i<n;i+=j) _mm512_store_pd((double *)&x[i],_mm512_load_pd((double *)&v.x[i]));
        } else if constexpr(std::is_same_v<float,Type> || std::is_same_v<complex<float>,Type>) {
            for(i=0;i<n;i+=j) _mm512_store_ps((float *)&x[i],_mm512_load_ps((float *)&v.x[i]));
        }
#endif
    }

    vector(vector<Type>&& v) {                                  // move constructor. need to modify v, no const
        //printf("vector construct &&\n");                      // class a = std::move(b);
        if(this == &v) return;
        n = v.n;
        if(x != NULL) free(x);
        x = v.x;
        v.n = 0;
        v.x = NULL;
    } 

    ~vector() {
        //printf("vector ~\n");
        if(x != NULL) { free(x); x = NULL; }
    }

    inline Type& operator [](int m) {
        return x[m];
    }

    inline Type operator [](int m) const {                // const here means not allow to change data of object
        return x[m];
    }

    //inline vector operator [](int m) const {
    //    return x[m];
    //}  

    vector<Type>& operator=(vector<Type>&& v) {                 // move assignment. need to modify v, no const
        //printf("vector= &&\n");                               // v3 = v1 - v2;
        if(this == &v) return *this;
        n = v.n;
        x = v.x;
        v.n = 0;
        v.x = NULL;
        return *this;
    }

    vector<Type>& operator=(const vector<Type>& v) {           // copy assignment
        int i,j;                                               // class a,b; a = b;
        //printf("vector= &\n");
        if(this == &v) return *this;
        n = v.n;
        if(x != NULL) { free(x); x = NULL; }
        x = allocate1D<Type>(n);
#if !defined(AVX) || AVX == 0
        memcpy(x,v.x,n*sizeof(Type));
#elif AVX512F == 0
        j = 32/sizeof(Type);
        if constexpr(std::is_same_v<double,Type> || std::is_same_v<complex<double>,Type>) {
            for(i=0;i<n;i+=j) _mm256_store_pd((double *)&x[i],_mm256_load_pd((double *)&v.x[i]));
        } else if constexpr(std::is_same_v<float,Type> || std::is_same_v<complex<float>,Type>) {
            for(i=0;i<n;i+=j) _mm256_store_ps((float *)&x[i],_mm256_load_ps((float *)&v.x[i]));
        }
#else
        j = 64/sizeof(Type);
        if constexpr(std::is_same_v<double,Type> || std::is_same_v<complex<double>,Type>) {
            for(i=0;i<n;i+=j) _mm512_store_pd((double *)&x[i],_mm512_load_pd((double *)&v.x[i]));
        } else if constexpr(std::is_same_v<float,Type> || std::is_same_v<complex<float>,Type>) {
            for(i=0;i<n;i+=j) _mm512_store_ps((float *)&x[i],_mm512_load_ps((float *)&v.x[i]));
        }
#endif
        return *this;
    }

    template <typename Type2>
    friend vector<Type2> operator+(const vector<Type2>&,const vector<Type2>&);
    template <typename Type2>
    friend void operator+=(const vector<Type2>&,const vector<Type2>&);
    template <typename Type2>
    friend vector<Type2> operator-(const vector<Type2>&,const vector<Type2>&);
    template <typename Type2>
    friend void operator-=(const vector<Type2>&,const vector<Type2>&);
    template <typename Type2>
    friend void swap(vector<Type2>&,vector<Type2>&);           // pass by reference
    template <typename Type2>
    friend Type2 dot(const vector<Type2>&,const vector<Type2>&);
    template <typename Type2>
    friend Type2 dot2(const vector<Type2>&,const vector<Type2>&);
    template <typename Type2>
    friend Type2 norm(const vector<Type2>&);
    template <typename Type2>
    friend ostream& operator<< (ostream&,const vector<Type2>&);
    template <typename Type2>
    friend vector<Type2> operator*(const Type2&,const vector<Type2>&);
    template <typename Type2>
    friend vector<Type2> operator*(const vector<Type2>&,const Type2&);
    template <typename Type2>
    friend void operator*=(const vector<Type2>&,const Type2&);
    //inline Type read(const int& m) const {
    //    return x[m];
    //}
    //inline void write(const int& m,const Type& a) {
    //    x[m] = a;
    //}
    friend class matrix<Type>;
};

template<typename Type2>
vector<Type2> operator+(const vector<Type2>& v1,const vector<Type2>& v2) {
    vector<Type2> temp(v1.n);
    int i,j;

#if !defined(AVX) || AVX == 0
    for(i=0;i<v1.n;i++) temp[i]=v1[i]+v2[i];
#elif AVX512F == 0
    j = 32/sizeof(Type2);
    if constexpr(std::is_same_v<double,Type2> || std::is_same_v<complex<double>,Type2>) {
        for(i=0;i<v1.n;i+=j) 
            _mm256_store_pd((double *)&temp.x[i],_mm256_add_pd(_mm256_load_pd((double *)&v1.x[i]),_mm256_load_pd((double *)&v2.x[i])));
    } else if constexpr(std::is_same_v<float,Type2> || std::is_same_v<complex<float>,Type2>) {
        for(i=0;i<v1.n;i+=j)
            _mm256_store_ps((float *)&temp.x[i],_mm256_add_ps(_mm256_load_ps((float *)&v1.x[i]),_mm256_load_ps((float *)&v2.x[i])));
    }
#else
    j = 64/sizeof(Type2);
    if constexpr(std::is_same_v<double,Type2> || std::is_same_v<complex<double>,Type2>) {
        for(i=0;i<v1.n;i+=j)
            _mm512_store_pd((double *)&temp.x[i],_mm512_add_pd(_mm512_load_pd((double *)&v1.x[i]),_mm512_load_pd((double *)&v2.x[i])));
    } else if constexpr(std::is_same_v<float,Type2> || std::is_same_v<complex<float>,Type2>) {
        for(i=0;i<v1.n;i+=j)
            _mm512_store_ps((float *)&temp.x[i],_mm512_add_ps(_mm512_load_ps((float *)&v1.x[i]),_mm512_load_ps((float *)&v2.x[i])));
    }
#endif
    return temp;
}

template<typename Type2>
void operator+=(const vector<Type2>& v1,const vector<Type2>& v2) {
    int i,j;
#if !defined(AVX) || AVX == 0
    for(i=0;i<v1.n;i++) v1[i]+=v2[i];
#elif AVX512F == 0
    j = 32/sizeof(Type2);
    if constexpr(std::is_same_v<double,Type2> || std::is_same_v<complex<double>,Type2>) {
        for(i=0;i<v1.n;i+=j) 
            _mm256_store_pd((double *)&v1.x[i],_mm256_add_pd(_mm256_load_pd((double *)&v1.x[i]),_mm256_load_pd((double *)&v2.x[i])));
    } else if constexpr(std::is_same_v<float,Type2> || std::is_same_v<complex<float>,Type2>) {
        for(i=0;i<v1.n;i+=j)
            _mm256_store_ps((float *)&v1.x[i],_mm256_add_ps(_mm256_load_ps((float *)&v1.x[i]),_mm256_load_ps((float *)&v2.x[i])));
    }
#else
    j = 64/sizeof(Type2);
    if constexpr(std::is_same_v<double,Type2> || std::is_same_v<complex<double>,Type2>) {
        for(i=0;i<v1.n;i+=j)
            _mm512_store_pd((double *)&v1.x[i],_mm512_add_pd(_mm512_load_pd((double *)&v1.x[i]),_mm512_load_pd((double *)&v2.x[i])));
    } else if constexpr(std::is_same_v<float,Type2> || std::is_same_v<complex<float>,Type2>) {
        for(i=0;i<v1.n;i+=j)
            _mm512_store_ps((float *)&v1.x[i],_mm512_add_ps(_mm512_load_ps((float *)&v1.x[i]),_mm512_load_ps((float *)&v2.x[i])));
    }
#endif
}

template<typename Type2>
vector<Type2> operator-(const vector<Type2>& v1,const vector<Type2>& v2) {
    vector<Type2> temp(v1.n);
    int i,j;

#if !defined(AVX) || AVX == 0
    for(i=0;i<v1.n;i++) temp[i]=v1[i]-v2[i];
#elif AVX512F == 0
    j = 32/sizeof(Type2);
    if constexpr(std::is_same_v<double,Type2> || std::is_same_v<complex<double>,Type2>) {
        for(i=0;i<v1.n;i+=j) 
            _mm256_store_pd((double *)&temp.x[i],_mm256_sub_pd(_mm256_load_pd((double *)&v1.x[i]),_mm256_load_pd((double *)&v2.x[i])));
    } else if constexpr(std::is_same_v<float,Type2> || std::is_same_v<complex<float>,Type2>) {
        for(i=0;i<v1.n;i+=j)
            _mm256_store_ps((float *)&temp.x[i],_mm256_sub_ps(_mm256_load_ps((float *)&v1.x[i]),_mm256_load_ps((float *)&v2.x[i])));
    }
#else
    j = 64/sizeof(Type2);
    if constexpr(std::is_same_v<double,Type2> || std::is_same_v<complex<double>,Type2>) {
        for(i=0;i<v1.n;i+=j)
            _mm512_store_pd((double *)&temp.x[i],_mm512_sub_pd(_mm512_load_pd((double *)&v1.x[i]),_mm512_load_pd((double *)&v2.x[i])));
    } else if constexpr(std::is_same_v<float,Type2> || std::is_same_v<complex<float>,Type2>) {
        for(i=0;i<v1.n;i+=j)
            _mm512_store_ps((float *)&temp.x[i],_mm512_sub_ps(_mm512_load_ps((float *)&v1.x[i]),_mm512_load_ps((float *)&v2.x[i])));
    }
#endif
    return temp;
}

template<typename Type2>
void operator-=(const vector<Type2>& v1,const vector<Type2>& v2) {
    int i,j;
#if !defined(AVX) || AVX == 0
    for(i=0;i<v1.n;i++) v1[i]-=v2[i];
#elif AVX512F == 0
    j = 32/sizeof(Type2);
    if constexpr(std::is_same_v<double,Type2> || std::is_same_v<complex<double>,Type2>) {
        for(i=0;i<v1.n;i+=j)
            _mm256_store_pd((double *)&v1.x[i],_mm256_sub_pd(_mm256_load_pd((double *)&v1.x[i]),_mm256_load_pd((double *)&v2.x[i])));
    } else if constexpr(std::is_same_v<float,Type2> || std::is_same_v<complex<float>,Type2>) {
        for(i=0;i<v1.n;i+=j)
            _mm256_store_ps((float *)&v1.x[i],_mm256_sub_ps(_mm256_load_ps((float *)&v1.x[i]),_mm256_load_ps((float *)&v2.x[i])));
    }
#else
    j = 64/sizeof(Type2);
    if constexpr(std::is_same_v<double,Type2> || std::is_same_v<complex<double>,Type2>) {
        for(i=0;i<v1.n;i+=j)
            _mm512_store_pd((double *)&v1.x[i],_mm512_sub_pd(_mm512_load_pd((double *)&v1.x[i]),_mm512_load_pd((double *)&v2.x[i])));
    } else if constexpr(std::is_same_v<float,Type2> || std::is_same_v<complex<float>,Type2>) {
        for(i=0;i<v1.n;i+=j)
            _mm512_store_ps((float *)&v1.x[i],_mm512_sub_ps(_mm512_load_ps((float *)&v1.x[i]),_mm512_load_ps((float *)&v2.x[i])));
    }
#endif
}

template<typename Type2>
void swap(vector<Type2>& v1,vector<Type2>& v2) {
    Type2 *temp;
    temp = v1.x;
    v1.x = v2.x;
    v2.x = temp;
    int m;
    m = v1.n;
    v1.n = v2.n;
    v2.n = m;
}

//  v1*v2.conjugate
template<typename Type2>
Type2 dot(const vector<Type2>& v1,const vector<Type2>& v2) {
    int i;
    Type2 tempsum(0.);
    if constexpr(std::is_same_v<float,Type2>||std::is_same_v<double,Type2>) {
        Type2 temp[v1.n];
        #pragma omp parallel for default(shared)
        for(i=0;i<v1.n;i++) temp[i]=v1[i]*v2[i];
        #pragma omp parallel for reduction(+:tempsum)
        for(i=0;i<v1.n;i++) tempsum+=temp[i];
    } else if constexpr(std::is_same_v<complex<float>,Type2>) {
        float temp_r[v1.n],temp_i[v1.n];
        float temp_r_sum = 0.;
        float temp_i_sum = 0.;
        #pragma omp parallel for default(shared)
        for(i=0;i<v1.n;i++) temp_r[i]=real(v1[i])*real(v2[i])+imag(v1[i])*imag(v2[i]);
        #pragma omp parallel for reduction(+:temp_r_sum)
        for(i=0;i<v1.n;i++) temp_r_sum+=temp_r[i];
        #pragma omp parallel for default(shared)
        for(i=0;i<v1.n;i++) temp_i[i]=imag(v1[i])*real(v2[i])-real(v1[i])*imag(v2[i]);
        #pragma omp parallel for reduction(+:temp_i_sum)
        for(i=0;i<v1.n;i++) temp_i_sum+=temp_i[i];
        tempsum = complex<float>(temp_r_sum,temp_i_sum);
    } else if constexpr(std::is_same_v<complex<double>,Type2>) {
        double temp_r[v1.n],temp_i[v1.n];
        double temp_r_sum = 0.;
        double temp_i_sum = 0.;
        #pragma omp parallel for default(shared)
        for(i=0;i<v1.n;i++) temp_r[i]=real(v1[i])*real(v2[i])+imag(v1[i])*imag(v2[i]);
        #pragma omp parallel for reduction(+:temp_r_sum)
        for(i=0;i<v1.n;i++) temp_r_sum+=temp_r[i];
        #pragma omp parallel for default(shared)
        for(i=0;i<v1.n;i++) temp_i[i]=imag(v1[i])*real(v2[i])-real(v1[i])*imag(v2[i]);
        #pragma omp parallel for reduction(+:temp_i_sum)
        for(i=0;i<v1.n;i++) temp_i_sum+=temp_i[i];
        tempsum = complex<double>(temp_r_sum,temp_i_sum);
    }
    return tempsum;
}

// v1*v2
template<typename Type2>
Type2 dot2(const vector<Type2>& v1,const vector<Type2>& v2) {
    int i;
    Type2 tempsum(0.);
    if constexpr(std::is_same_v<float,Type2>||std::is_same_v<double,Type2>) {
        Type2 temp[v1.n];
        #pragma omp parallel for default(shared)
        for(i=0;i<v1.n;i++) temp[i]=v1[i]*v2[i];
        #pragma omp parallel for reduction(+:tempsum)
        for(i=0;i<v1.n;i++) tempsum+=temp[i];
    } else if constexpr(std::is_same_v<complex<float>,Type2>) {
        float temp_r[v1.n],temp_i[v1.n];
        float temp_r_sum = 0.;
        float temp_i_sum = 0.;
        #pragma omp parallel for default(shared)
        for(i=0;i<v1.n;i++) temp_r[i]=real(v1[i])*real(v2[i])-imag(v1[i])*imag(v2[i]);
        #pragma omp parallel for reduction(+:temp_r_sum)
        for(i=0;i<v1.n;i++) temp_r_sum+=temp_r[i];
        #pragma omp parallel for default(shared)
        for(i=0;i<v1.n;i++) temp_i[i]=imag(v1[i])*real(v2[i])+real(v1[i])*imag(v2[i]);
        #pragma omp parallel for reduction(+:temp_i_sum)
        for(i=0;i<v1.n;i++) temp_i_sum+=temp_i[i];
        tempsum = complex<float>(temp_r_sum,temp_i_sum);
    } else if constexpr(std::is_same_v<complex<double>,Type2>) {
        double temp_r[v1.n],temp_i[v1.n];
        double temp_r_sum = 0.;
        double temp_i_sum = 0.;
        #pragma omp parallel for default(shared)
        for(i=0;i<v1.n;i++) temp_r[i]=real(v1[i])*real(v2[i])-imag(v1[i])*imag(v2[i]);
        #pragma omp parallel for reduction(+:temp_r_sum)
        for(i=0;i<v1.n;i++) temp_r_sum+=temp_r[i];
        #pragma omp parallel for default(shared)
        for(i=0;i<v1.n;i++) temp_i[i]=imag(v1[i])*real(v2[i])+real(v1[i])*imag(v2[i]);
        #pragma omp parallel for reduction(+:temp_i_sum)
        for(i=0;i<v1.n;i++) temp_i_sum+=temp_i[i];
        tempsum = complex<double>(temp_r_sum,temp_i_sum);
    }
    return tempsum;
}

template<typename Type2> 
Type2 norm(const vector<Type2>& v) {
    int i;
    Type2 val;
    if constexpr(std::is_same_v<float,Type2> || std::is_same_v<double,Type2>) {
        Type2 temp[v.n],tempsum = 0.;
        #pragma omp parallel for default(shared)
        for(i=0;i<v.n;i++) temp[i]=v[i]*v[i];
        #pragma omp parallel for reduction(+:tempsum)
        for(i=0;i<v.n;i++) tempsum+=temp[i];
        tempsum = sqrt(tempsum);
        val = tempsum;
    } else if constexpr(std::is_same_v<complex<float>,Type2>) {
        float temp[v.n],tempsum = 0.;
        #pragma omp parallel for default(shared)
        for(i=0;i<v.n;i++) temp[i]=norm(v[i]);
        #pragma omp parallel for reduction(+:tempsum)
        for(i=0;i<v.n;i++) tempsum+=temp[i];
        tempsum = sqrt(tempsum);
        val = tempsum; 
    } else if constexpr(std::is_same_v<complex<double>,Type2>) {
        double temp[v.n],tempsum = 0.;
        #pragma omp parallel for default(shared)
        for(i=0;i<v.n;i++) temp[i]=norm(v[i]);
        #pragma omp parallel for reduction(+:tempsum)
        for(i=0;i<v.n;i++) tempsum+=temp[i];
        tempsum = sqrt(tempsum);
        val = tempsum;
    }
    return val;
}

template<typename Type2>
ostream& operator<< (ostream& s,const vector<Type2>& v) {
    for(int i=0;i<v.n;i++) s << v[i] << "      ";
    s << endl;
    return s;
}

template<typename Type2>
vector<Type2> operator*(const Type2& a,const vector<Type2>& v) {
    vector<Type2> temp(v.n);
    int i,j;

#if AVX512F > 0
    typedef typename std::conditional<std::is_same_v<double,Type2> || std::is_same_v<complex<double>,Type2>,__m512d,__m512>::type avxtype;
    alignas(ALIGN) avxtype mw01_vecmat_a;
#elif AVX > 0
    typedef typename std::conditional<std::is_same_v<double,Type2> || std::is_same_v<complex<double>,Type2>,__m256d,__m256>::type avxtype;
    alignas(ALIGN) avxtype mw01_vecmat_a;
#endif

#if !defined(AVX) || AVX == 0
    for(i=0;i<v.n;i++) temp[i] = v[i]*a;
#elif AVX512F == 0
    j = 32/sizeof(Type2);
    if constexpr(std::is_same_v<double,Type2>) {
        mw01_vecmat_a = _mm256_set1_pd(a);
        for(i=0;i<v.n;i+=j)
            _mm256_store_pd((double *)&temp.x[i],_mm256_mul_pd(mw01_vecmat_a,_mm256_load_pd((double *)&v.x[i]))); 
    } else if constexpr(std::is_same_v<complex<double>,Type2>) {
        mw01_vecmat_a = _mm256_setr_pd(real(a),imag(a),real(a),imag(a));
        for(i=0;i<v.n;i+=j)
            _mm256_store_pd((double *)&temp.x[i],complex_mul_256register(mw01_vecmat_a,_mm256_load_pd((double *)&v.x[i])));
    } else if constexpr(std::is_same_v<float,Type2>) {
        mw01_vecmat_a = _mm256_set1_ps(a);
        for(i=0;i<v.n;i+=j)
            _mm256_store_ps((float *)&temp.x[i],_mm256_mul_ps(mw01_vecmat_a,_mm256_load_ps((float *)&v.x[i])));
    } else if constexpr(std::is_same_v<complex<float>,Type2>) {
        mw01_vecmat_a = _mm256_setr_ps(real(a),imag(a),real(a),imag(a),real(a),imag(a),real(a),imag(a));
        for(i=0;i<v.n;i+=j)
            _mm256_store_ps((float *)&temp.x[i],complex_mul_256register(mw01_vecmat_a,_mm256_load_ps((float *)&v.x[i])));
    }
#else
    j = 64/sizeof(Type2);
    if constexpr(std::is_same_v<double,Type2>) {
        mw01_vecmat_a = _mm512_set1_pd(a);
        for(i=0;i<v.n;i+=j)
            _mm512_store_pd((double *)&temp.x[i],_mm512_mul_pd(mw01_vecmat_a,_mm512_load_pd((double *)&v.x[i])));
    } else if constexpr(std::is_same_v<complex<double>,Type2>) {
        mw01_vecmat_a = _mm512_setr_pd(real(a),imag(a),real(a),imag(a),real(a),imag(a),real(a),imag(a));
        for(i=0;i<v.n;i+=j)
            _mm512_store_pd((double *)&temp.x[i],complex_mul_512register(mw01_vecmat_a,_mm512_load_pd((double *)&v.x[i])));
    } else if constexpr(std::is_same_v<float,Type2>) {
        mw01_vecmat_a = _mm512_set1_ps(a);
        for(i=0;i<v.n;i+=j)
            _mm512_store_ps((float *)&temp.x[i],_mm512_mul_ps(mw01_vecmat_a,_mm512_load_ps((float *)&v.x[i])));
    } else if constexpr(std::is_same_v<complex<float>,Type2>) {
        mw01_vecmat_a = _mm512_setr_ps(real(a),imag(a),real(a),imag(a),real(a),imag(a),real(a),imag(a),real(a),imag(a),real(a),imag(a),real(a),imag(a),real(a),imag(a));
        for(i=0;i<v.n;i+=j)
            _mm512_store_ps((float *)&temp.x[i],complex_mul_512register(mw01_vecmat_a,_mm512_load_ps((float *)&v.x[i])));
    }
#endif
    return temp;
}

template<typename Type2>
vector<Type2> operator*(const vector<Type2>& v,const Type2& a) {
    return a*v;
}

template<typename Type2>
void operator*=(const vector<Type2>& v,const Type2& a) {
    int i,j;

#if AVX512F > 0
    typedef typename std::conditional<std::is_same_v<double,Type2> || std::is_same_v<complex<double>,Type2>,__m512d,__m512>::type avxtype;
    alignas(ALIGN) avxtype mw01_vecmat_a;
#elif AVX > 0
    typedef typename std::conditional<std::is_same_v<double,Type2> || std::is_same_v<complex<double>,Type2>,__m256d,__m256>::type avxtype;
    alignas(ALIGN) avxtype mw01_vecmat_a;
#endif

#if !defined(AVX) || AVX == 0
    for(i=0;i<v.n;i++) v[i] = v[i]*a;
#elif AVX512F == 0
    j = 32/sizeof(Type2);
    if constexpr(std::is_same_v<double,Type2>) {
        mw01_vecmat_a = _mm256_set1_pd(a);
        for(i=0;i<v.n;i+=j)
            _mm256_store_pd((double *)&v.x[i],_mm256_mul_pd(mw01_vecmat_a,_mm256_load_pd((double *)&v.x[i]))); 
    } else if constexpr(std::is_same_v<complex<double>,Type2>) {
        mw01_vecmat_a = _mm256_setr_pd(real(a),imag(a),real(a),imag(a));
        for(i=0;i<v.n;i+=j)
            _mm256_store_pd((double *)&v.x[i],complex_mul_256register(mw01_vecmat_a,_mm256_load_pd((double *)&v.x[i])));
    } else if constexpr(std::is_same_v<float,Type2>) {
        mw01_vecmat_a = _mm256_set1_ps(a);
        for(i=0;i<v.n;i+=j)
            _mm256_store_ps((float *)&v.x[i],_mm256_mul_ps(mw01_vecmat_a,_mm256_load_ps((float *)&v.x[i])));
    } else if constexpr(std::is_same_v<complex<float>,Type2>) {
        mw01_vecmat_a = _mm256_setr_ps(real(a),imag(a),real(a),imag(a),real(a),imag(a),real(a),imag(a));
        for(i=0;i<v.n;i+=j)
            _mm256_store_ps((float *)&v.x[i],complex_mul_256register(mw01_vecmat_a,_mm256_load_ps((float *)&v.x[i])));
    }
#else
    j = 64/sizeof(Type2);
    if constexpr(std::is_same_v<double,Type2>) {
        mw01_vecmat_a = _mm512_set1_pd(a);
        for(i=0;i<v.n;i+=j)
            _mm512_store_pd((double *)&v.x[i],_mm512_mul_pd(mw01_vecmat_a,_mm512_load_pd((double *)&v.x[i])));
    } else if constexpr(std::is_same_v<complex<double>,Type2>) {
        mw01_vecmat_a = _mm512_setr_pd(real(a),imag(a),real(a),imag(a),real(a),imag(a),real(a),imag(a));
        for(i=0;i<v.n;i+=j)
            _mm512_store_pd((double *)&v.x[i],complex_mul_512register(mw01_vecmat_a,_mm512_load_pd((double *)&v.x[i])));
    } else if constexpr(std::is_same_v<float,Type2>) {
        mw01_vecmat_a = _mm512_set1_ps(a);
        for(i=0;i<v.n;i+=j)
            _mm512_store_ps((float *)&v.x[i],_mm512_mul_ps(mw01_vecmat_a,_mm512_load_ps((float *)&v.x[i])));
    } else if constexpr(std::is_same_v<complex<float>,Type2>) {
        mw01_vecmat_a = _mm512_setr_ps(real(a),imag(a),real(a),imag(a),real(a),imag(a),real(a),imag(a),real(a),imag(a),real(a),imag(a),real(a),imag(a),real(a),imag(a));
        for(i=0;i<v.n;i+=j)
            _mm512_store_ps((float *)&v.x[i],complex_mul_512register(mw01_vecmat_a,_mm512_load_ps((float *)&v.x[i])));
    }
#endif
}




template <class Type>
class matrix {
    int m = 0;
    int n = 0;
    vector<Type>* y = NULL;
  public:
    matrix() {
        for(int k=0;k<m;k++) if(y[k].x != NULL) { free(y[k].x); y[k].x = NULL; }
        if(y != NULL) { delete[] y; y = NULL; }
        m=0;
        n=0;
    }
    matrix(int i,int j) {
        for(int k=0;k<m;k++) if(y[k].x != NULL) { free(y[k].x); y[k].x = NULL; }
        if(y != NULL) { delete[] y; y = NULL; }
        m=i;
        n=j;
        y = new vector<Type>[m];                            
        for(int k=0;k<m;k++) y[k] = vector<Type>(n);        
    }
    matrix(int i,int j,Type* a) {
        for(int k=0;k<m;k++) if(y[k].x != NULL) { free(y[k].x); y[k].x = NULL; }
        if(y != NULL) { delete[] y; y = NULL; }
        m=i;
        n=j;
        y = new vector<Type>[m];
        for(int k=0;k<m;k++,a+=n) y[k] = vector<Type>(n,a);
    }

    matrix(const matrix<Type>& M) {                       // copy constructor                      
        printf("matrix &\n");
        if(this == &M) return;
        for(int k=0;k<m;k++) if(y[k].x != NULL) { free(y[k].x); y[k].x = NULL; }
        if(y != NULL) { delete[] y; y = NULL; }
        m = M.m;
        n = M.n;
        y = new vector<Type>[m];
        for(int k=0;k<m;k++) y[k] = M[k];
    }

    matrix(matrix<Type>&& M) {                            // move constructor
        printf("matrix move construct\n");
        if(this == &M) return;
        for(int k=0;k<m;k++) if(y[k].x != NULL) { free(y[k].x); y[k].x = NULL; }
        if(y != NULL) { delete[] y; y = NULL; }
        m = M.m;
        n = M.n;
        y = M.y;
        M.m = 0;
        M.n = 0;
        M.y = NULL;
    }

    ~matrix() {
        printf("matrix ~\n");
        for(int k=0;k<m;k++) if(y[k].x != NULL) { free(y[k].x); y[k].x = NULL; }
        if(y != NULL) { delete[] y; y = NULL; }
    }

    inline vector<Type>& operator [](int m) {
        return y[m];
    } 

    inline vector<Type> operator [](int m) const {
       return y[m];
    }

    matrix<Type>& operator=(matrix<Type>&& M) {           // move assignment
        printf("matrix move assign\n");
        if(this == &M) return *this;
        for(int k=0;k<m;k++) if(y[k].x != NULL) { free(y[k].x); y[k].x = NULL; }
        if(y != NULL) { delete[] y; y = NULL; }
        m = M.m;
        n = M.n;
        y = M.y;
        M.m = 0;
        M.n = 0;
        M.y = NULL;
        return *this;
    }

    matrix<Type>& operator=(const matrix<Type>& M) {                    //  copy assignment
        printf("matrix =\n");
        if(this == &M) return *this;
        for(int k=0;k<m;k++) if(y[k].x != NULL) { free(y[k].x); y[k].x = NULL; }
        if(y != NULL) { delete[] y; y = NULL; }
        m = M.m;
        n = M.n;
        y = new vector<Type>[m];
        for(int k=0;k<m;k++) y[k] = M[k];
        return *this;
    }

    template<typename Type2>
    friend matrix<Type2> transpose(const matrix<Type2>&);                      
    template<typename Type2>
    friend matrix<Type2> conjugate_transpose(const matrix<Type2>&);            
    template<typename Type2>
    friend matrix<Type2> operator+(const matrix<Type2>&,const matrix<Type2>&);
    template<typename Type2>
    friend void operator+=(const matrix<Type2>&,const matrix<Type2>&);
    template<typename Type2>
    friend matrix<Type2> operator-(const matrix<Type2>&,const matrix<Type2>&);
    template<typename Type2>
    friend void operator-=(const matrix<Type2>&,const matrix<Type2>&);
    template<typename Type2>
    friend matrix<Type2> operator*(const matrix<Type2>&,const matrix<Type2>&);
    template<typename Type2>
    friend void operator*=(const matrix<Type2>&,const matrix<Type2>&);
    template<typename Type2>
    friend matrix<Type2> LU(matrix<Type2>,int *,int);
    template<typename Type2>
    friend matrix<Type2> QR(matrix<Type2>,int *);
    template<typename Type2>
    friend matrix<Type2> CHOLESKY(matrix<Type2>,int *,int);
    template<typename Type2>
    friend matrix<Type2> inverse_LowerTriangular(const matrix<Type2>&);
    template<typename Type2>
    friend matrix<Type2> inverse_LU(const matrix<Type2>&,int *);
    template<typename Type2>
    friend matrix<Type2> inverse_QR(const matrix<Type2>&,int *);
    template<typename Type2>
    friend matrix<Type2> inverse_CHOLESKY(const matrix<Type2>&,int *,int);
    template<typename Type2>
    friend matrix<Type2> inverse_GAUSS(const matrix<Type2>&,int);
    template<typename Type2>
    friend vector<Type2> solveAxb_LU(const matrix<Type2>&,int *,const vector<Type2>&);
    template<typename Type2>
    friend vector<Type2> solveAxb_GAUSS(const matrix<Type2>&,const vector<Type2>&,int);
    template<typename Type2>
    friend vector<Type2> solveAxb_CHOLESKY(const matrix<Type2>&,int *,const matrix<Type2>&,const vector<Type2>&);
    template<typename Type2>
    friend matrix<Type2> operator*(const Type2&,const matrix<Type2>&);
    template<typename Type2>
    friend matrix<Type2> operator*(const matrix<Type2>&,const Type2&);
    template<typename Type2>
    friend void operator*=(const matrix<Type2>&,const Type2&);
    template<typename Type2>
    friend vector<Type2> operator*(const matrix<Type2>&,const vector<Type2>&);
    template<typename Type2>
    friend ostream &operator<<(ostream& s,matrix<Type2>& M);
    inline vector<Type> readvector(const int& m) {
        return y[m]; 
    }
    inline void writevector(const int& m,const vector<Type>& v) {
        y[m] = v;
    }
    inline vector<Type> readcolvector(const int& n) const {
        int i;
        vector<Type> v(m);
        for(i=0;i<m;i++) v[i] = y[i][n];
        return v;
    }
    inline void writecolvector(const int& n,const vector<Type>& v) {
        int i;
        for(i=0;i<m;i++) y[i][n] = v[i];
    }
    void *address(const int& m) {
        return y[m].x; 
    }
    friend class vector<Type>;
};

template<typename Type2>
matrix<Type2> transpose(const matrix<Type2>& M) {
    int i,j;
    matrix<Type2> temp(M.n,M.m);
    #pragma omp parallel for default(shared) private(j)
    for(i=0;i<M.n;i++)
    for(j=0;j<M.m;j++)
        temp[i][j] = M[j][i];
    return temp;
}

template<typename Type2>
matrix<Type2> conjugate_transpose(const matrix<Type2>& M) {
    int i,j;
    matrix<Type2> temp(M.n,M.m);
    #pragma omp parallel for default(shared) private(j)
    for(i=0;i<M.n;i++)
    for(j=0;j<M.m;j++) {
        if constexpr(std::is_same_v<float,Type2> || std::is_same_v<double,Type2>) {
            temp[i][j] = M[j][i]; 
        } else if constexpr(std::is_same_v<complex<float>,Type2> || std::is_same_v<complex<double>,Type2>) {
            temp[i][j] = conj(M[j][i]);
        }
    }
    return temp;
}

template<typename Type2>
matrix<Type2> operator+(const matrix<Type2>& M1,const matrix<Type2>& M2) {
    matrix<Type2> temp(M1.m,M1.n);
    for(int i=0;i<M1.m;i++) temp[i] = M1[i] + M2[i];
    return temp;
}

template<typename Type2>
void operator+=(const matrix<Type2>& M1,const matrix<Type2>& M2) {
    for(int i=0;i<M1.m;i++) M1[i] += M2[i];
}

template<typename Type2>
matrix<Type2> operator-(const matrix<Type2>& M1,const matrix<Type2>& M2) {
    matrix<Type2> temp(M1.m,M1.n);
    for(int i=0;i<M1.m;i++) temp[i] = M1[i] - M2[i];
    return temp;
}

template<typename Type2>
void operator-=(const matrix<Type2>& M1,const matrix<Type2>& M2) {
    for(int i=0;i<M1.m;i++) M1[i] -= M2[i];
}

template<typename Type2>
matrix<Type2> operator*(const matrix<Type2>& M1,const matrix<Type2>& M2) {
    matrix<Type2> temp(M1.m,M2.n);
    matrix<Type2> Mtran = transpose(M2);
    for(int i=0;i<M1.m;i++)
    for(int j=0;j<M2.n;j++)
        temp[i][j] = dot2(M1[i],Mtran[j]);
    return temp;
}

template<typename Type2>
void operator*=(const matrix<Type2>& M1,const matrix<Type2>& M2) {
    matrix<Type2> temp(M1.m,M2.n);
    matrix<Type2> Mtran = transpose(M2);
    for(int i=0;i<M1.m;i++)
    for(int j=0;j<M2.n;j++)
        temp[i][j] = dot2(M1[i],Mtran[j]);
    M1 = temp;                                                // possible error here?
}

template<typename Type2>
matrix<Type2> LU(matrix<Type2> A,int* P,int pivot) {    // should not modify the calling M, so not by reference
    int i,j,k,l;
    double a[A.m];
    double maxabs;
    //matrix<Type2> A = M; 

    #pragma omp parallel for default(shared)
    for(j=0;j<=A.m;j++) P[j] = j;

    for(j=0;j<A.m;j++) {
        if(pivot == 0) {
            for(k=j;k<A.m;k++) if(A[k][j] != Type2(0.)) break;
        } else if(pivot == 1) {
            maxabs = 0.;
            #pragma omp parallel for default(shared)
            for(l=j;l<A.m;l++) a[l] = fabs(A[l][j]);
            for(l=j;l<A.m;l++) if(a[l] > maxabs) { maxabs = a[l]; k = l; }
        }
        if(k != j) {
            //pivoting P
            l = P[j];
            P[j] = P[k];
            P[k] = l;
            //pivoting rows of A
            swap(A.y[j],A.y[k]);
            //counting pivots starting from N (for determinant)
            P[A.m]++;
        }
        for(k=j+1;k<A.m;k++) {
            A[k][j] = A[k][j]/A[j][j];
            #pragma omp parallel for default(shared)
            for(l=j+1;l<A.m;l++)
                A[k][l] = A[k][l] - A[k][j]*A[j][l];
        }
    }
    return A;
}

template<typename Type2>
matrix<Type2> QR(matrix<Type2> M,int* P) {    
    int i,j,k,l;
    int pivot = 1;
    matrix<Type2> A = M; 

    return A;
}

//   Hermitian positive definite A,   PT = P-1
//    P*A*PT = L*LH,    A = PT*L*LH*P = PT*L*(PT*L)^H
//   pivotted code:         
template<typename Type2>
matrix<Type2> CHOLESKY(matrix<Type2> Q,int* P,int pivot) {
    int i,j,k,l;
    matrix<Type2> L(Q.m,Q.m);
    Type2 sum[Q.m];
    Type2 tempsum;
    Type2 ZERO(0.);
    double a[Q.m];
    double maxabs;

    #pragma omp parallel for default(shared)
    for(j=0;j<Q.m;j++) P[j] = j;

 
    for(j=0;j<Q.m;j++) {
        if(pivot == 0) {
            k = j;
        } else if(pivot == 1) {                                        // partial pivoting
            maxabs = 0.;
            #pragma omp parallel for default(shared)
            for(l=j;l<Q.m;l++) a[l] = fabs(Q[l][l]);
            for(l=j;l<Q.m;l++) if(a[l] > maxabs) { maxabs = a[l]; k = l; }
        }
        if(k != j) {
            l = P[j];
            P[j] = P[k];
            P[k] = l;
            swap(Q.y[j],Q.y[k]);
            #pragma omp parallel for default(shared) private(tempsum)
            for(l=0;l<Q.m;l++) {
                tempsum = Q[l][j];
                Q[l][j] = Q[l][k];
                Q[l][k] = tempsum;
            }  
            P[Q.m]++;
        }
        for(i=0;i<=j;i++) {
            tempsum = 0.;
            #pragma omp parallel for default(shared)
            for(k=0;k<i;k++)
                if(L[j][k] == ZERO || L[i][k] == ZERO) {
                    sum[k] = 0.;
                } else {
                    if constexpr(std::is_same_v<float,Type2> || std::is_same_v<double,Type2>) {
                        sum[k] = L[j][k]*L[i][k];
                    } else if constexpr(std::is_same_v<complex<float>,Type2> || std::is_same_v<complex<double>,Type2>) {
                        sum[k] = L[j][k]*conj(L[i][k]);
                    }
                }
            for(k=0;k<i;k++) if(sum[k] != ZERO) tempsum += sum[k];
            if(j == i) {
                if constexpr(std::is_same_v<float,Type2> || std::is_same_v<double,Type2>) {
                    L[j][i] = sqrt(Q[j][j]-tempsum);
                } else if constexpr(std::is_same_v<complex<float>,Type2> || std::is_same_v<complex<double>,Type2>) {
                    L[j][i] = Type2(sqrt(real(Q[j][j])-real(tempsum)),0.);
                }
            } else {
                L[j][i] = (Q[j][i]-tempsum)/L[i][i]; 
            }
        }
    }
    return L;
}

template<typename Type2>
matrix<Type2> inverse_LowerTriangular(const matrix<Type2>& L) {
    int i,j,k;
    matrix<Type2> L_INV(L.m,L.m);

    for(i=0;i<L.m;i++) L_INV[i][i] = Type2(1.)/L[i][i];
    for(j=0;j<L.m;j++)
        for(i=j+1;i<L.m;i++) {
            L_INV[i][j] = Type2(0.);
            for(k=j;k<i;k++) L_INV[i][j] = L_INV[i][j]-L[i][k]*L_INV[k][j];
            L_INV[i][j] = L_INV[i][j]/L[i][i];
        }
    for(i=0;i<L.m;i++)
        for(j=i+1;j<L.m;j++) L_INV[i][j] = Type2(0.);
    return L_INV;
}

template<typename Type2>
matrix<Type2> inverse_LU(const matrix<Type2>& LU,int* P) {
    matrix<Type2> inv(LU.m,LU.m);
    for(int i=0;i<LU.m;i++) {
        for(int j=0;j<LU.m;j++) {
            if(P[j] == i) inv[j][i] = Type2(1.0); else inv[j][i] = Type2(0.0);
            #pragma omp parallel for default(shared)
            for(int k=0;k<j;k++)
                inv[j][i] = inv[j][i] - LU[j][k]*inv[k][i];
        }
        for(int j=LU.m-1;j>=0;j--) {
            #pragma omp parallel for default(shared)
            for(int k=j+1;k<LU.m;k++)
                inv[j][i] = inv[j][i]-LU[j][k]*inv[k][i];
            inv[j][i] = inv[j][i]/LU[j][j];
        }
    }
    return inv;
}

template<typename Type2>
matrix<Type2> inverse_QR(const matrix<Type2>& LU,int* P) {
    matrix<Type2> inv(LU.m,LU.m);
    return inv;
}

// A = PT*LLH*P = (PT*L)*(PT*L)^H
// A-1 = PT*(LLH)^-1*P = PT*(L^-1)^H*(L^-1)*P
// Ax = b, PT*LLH*Px = b, LLH*Px = Pb
template<typename Type2>
matrix<Type2> inverse_CHOLESKY(const matrix<Type2>& L,int* P,int flag) {    // flag = 0: input is L,  = 1: input is L-1
    matrix<Type2> L_INV(L.m,L.m);
    matrix<Type2> temp(L.m,L.m),temp2(L.m,L.m);
    int PT[L.m];
    int i,j,k;

    #pragma omp parallel for default(shared)
    for(j=0;j<L.m;j++) PT[P[j]] = j;

    if(flag == 0) {
        L_INV = inverse_LowerTriangular(L);
    } else {
        L_INV = L;
    }

    temp = conjugate_transpose(L_INV)*L_INV;

    for(i=0;i<L.m;i++)
    for(j=0;j<L.m;j++)
        temp2[P[i]][P[j]] = temp[i][j];

    return temp2; 
}

// By Gauss elimination
template<typename Type2>
matrix<Type2> inverse_GAUSS(const matrix<Type2>& M,int pivot) {
    int i,j,k,l;
    matrix<Type2> inv(M.m,M.m);
    matrix<Type2> temp(M.m,2*M.m);
    double a[M.m];
    double maxabs;

    #pragma omp parallel for default(shared) private(i)
    for(j=0;j<M.m;j++) {
        for(i=0;i<M.m;i++) temp[j][i] = M[j][i];
        temp[j][j+M.m] = Type2(1.);
    }
    for(j=0;j<M.m;j++) {
        if(pivot == 0) {
            for(k=j;k<M.m;k++) if(temp[k][j] != Type2(0.)) break;
        } else if(pivot == 1) {
            maxabs = 0.;
            #pragma omp parallel for default(shared)
            for(l=j;l<M.m;l++) a[l] = fabs(temp[l][j]);
            for(l=j;l<M.m;l++) if(a[l] > maxabs) { maxabs = a[l]; k = l; }
        }
        if(k == M.m) {
            cout << "Determinant = 0" << endl;
            return inv;
        }
        if(k != j) swap(temp.y[j],temp.y[k]);
        if(temp[j][j] != Type2(1.)) {
            #pragma omp parallel for default(shared)
            for(i=j+1;i<2*M.m;i++) temp[j][i] = temp[j][i]/temp[j][j];
            temp[j][j] = Type2(1.);
        }
        for(k=1;k<M.m;k++) {
            l=(j+k)%M.m;
            if(temp[l][j] != Type2(0.)) {
                #pragma omp parallel for default(shared)
                for(i=j+1;i<2*M.m;i++) temp[l][i] = temp[l][i]-temp[l][j]*temp[j][i];
                temp[l][j] = Type2(0.);
            }
        }
    }

    #pragma omp parallel for default(shared) private(i)
    for(j=0;j<M.m;j++)
    for(i=0;i<M.m;i++)
        inv[j][i] = temp[j][i+M.m];
    return inv;
}

// From wiki:   https://en.wikipedia.org/wiki/LU_decomposition
template<typename Type2>
vector<Type2> solveAxb_LU(const matrix<Type2>& LU,int* P,const vector<Type2>& b) {      
    vector<Type2> x(LU.m);

    for(int j=0;j<LU.m;j++) {
        x[j] = b[P[j]]; 
        #pragma omp parallel for default(shared)
        for(int k=0;k<j;k++)
            x[j] = x[j]-LU[j][k]*x[k];
    }
    for(int j=LU.m-1;j>=0;j--) {
        #pragma omp parallel for default(shared)
        for(int k=j+1;k<LU.m;k++)
            x[j] = x[j]-LU[j][k]*x[k];
        x[j] = x[j]/LU[j][j];
    }
    return x;
}

template<typename Type2>
vector<Type2> solveAxb_GAUSS(const matrix<Type2>& M,const vector<Type2>& b,int pivot) {
    int i,j,k,l,kk,ll;
    vector<Type2> x(M.m);
    matrix<Type2> temp(M.m,M.m+1);
    Type2 c;
    int permute[M.m];
    double a[M.m];
    double maxabs;

    for(l=0;l<M.m;l++) permute[l] = l;

    #pragma omp parallel for default(shared) private(i)
    for(j=0;j<M.m;j++) {                                                   
        for(i=0;i<M.m;i++) temp[j][i] = M[j][i];                   // j-th row; i-th col
        temp[j][M.m] = b[j];
    }
    printf("M.m = %d\n",M.m);

    for(j=0;j<M.m;j++) {
        if(pivot == 0) {
            for(k=j;k<M.m;k++) if(temp[k][j] != Type2(0.)) break;         // k-th row; j-th col
        } else if(pivot == 1) {
            maxabs = 0.;
            #pragma omp parallel for default(shared)
            for(l=j;l<M.m;l++) a[l] = fabs(temp[l][j]);               // l-th row; j-th col
            for(l=j;l<M.m;l++) if(a[l] > maxabs) { maxabs = a[l]; k = l; }
        } else if(pivot == 2) {
            maxabs = 0.;
            for(l=j;l<M.m;l++) {              // l-th row; ll-th col
                #pragma omp parallel for default(shared)
                for(ll=j;ll<M.m;ll++) a[ll] = fabs(temp[l][ll]);
                for(ll=j;ll<M.m;ll++) if(a[ll] > maxabs) { maxabs = a[ll]; k = l; kk = ll; }
            }

            // start swapping columns
            if(kk != j) {
                ll = permute[kk];
                permute[kk] = permute[j];
                permute[j] = ll;
                for(l=0;l<M.m;l++) {
                    c = temp[l][kk];
                    temp[l][kk] = temp[l][j];
                    temp[l][j] = c;
                }
            }
        }
        if(k == M.m) {
            cout << "Determinant = 0" << endl;
            return x;
        }
        if(k != j) swap(temp.y[j],temp.y[k]);

        if(temp[j][j] != Type2(1.)) {
            #pragma omp parallel for default(shared)
            for(i=j+1;i<M.m+1;i++) temp[j][i] = temp[j][i]/temp[j][j];
            temp[j][j] = Type2(1.);
        }
        for(k=1;k<M.m;k++) {
            l=(j+k)%M.m;
            if(temp[l][j] != Type2(0.)) {
                #pragma omp parallel for default(shared)
                for(i=j+1;i<M.m+1;i++) temp[l][i] = temp[l][i]-temp[l][j]*temp[j][i];
                temp[l][j] = Type2(0.);
            }
        }
    }

    #pragma omp parallel for default(shared)
    for(j=0;j<M.m;j++)
        x[permute[j]] = temp[j][M.m];
        //x.write(j,temp.read(j,M.m));
    return x;
}

// A = PT*LLH*P
// Ax = b, PT*LLH*Px = b, LLH*Px = Pb, [ Ly=Pb, (LH)Px = y ]
template<typename Type2>
vector<Type2> solveAxb_CHOLESKY(const matrix<Type2>& L,int* P,const matrix<Type2>& LH,const vector<Type2>& b) {
    vector<Type2> y(L.m),x(L.m),Px(L.m);
    int PT[L.m];
    int i,j,k;

    #pragma omp parallel for default(shared)
    for(j=0;j<L.m;j++) PT[P[j]] = j;

    for(j=0;j<L.m;j++) {
        y[j] = b[P[j]];
        for(i=0;i<j;i++) y[j] = y[j]-L[j][i]*y[i];
        y[j] = y[j]/L[j][j]; 
    }

    for(j=L.m-1;j>=0;j--) {
        Px[j] = y[j];
        for(i=j+1;i<L.m;i++) Px[j] = Px[j]-LH[j][i]*Px[i];
        Px[j] = Px[j]/LH[j][j];
    }
    #pragma omp parallel for default(shared)
    for(j=0;j<L.m;j++) x[j] = Px[PT[j]];

    return x;
}

template<typename Type2>
matrix<Type2> operator*(const Type2& a,const matrix<Type2>& M) {
    int i;
    matrix<Type2> temp(M.m,M.n);
    for(i=0;i<M.m;i++) temp[i] = a*M[i];
    return temp;
}

template<typename Type2>
matrix<Type2> operator*(const matrix<Type2>& M,const Type2& a) {
    return a*M;
}

template<typename Type2>
void operator*=(const matrix<Type2>& M,const Type2& a) {
    int i;
    matrix<Type2> temp(M.m,M.n);
    for(i=0;i<M.m;i++) M[i] *= a;
}

template<typename Type2>
vector<Type2> operator*(const matrix<Type2>& M,const vector<Type2>& v) {   // without &, use implicit copy/move ?
    vector<Type2> temp(M.m);
    for(int i=0;i<M.m;i++) temp[i] = dot2(M[i],v);
    return temp;
}

template<typename Type2>
ostream &operator<<(ostream& s,matrix<Type2>& M) {
    for(int k=0;k<M.m;k++) s << M[k];
    return s;
}








