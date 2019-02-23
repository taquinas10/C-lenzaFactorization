#include<math.h>
#include<stdio.h>
#include<gmp.h>
#include<stdarg.h>
#include<time.h>


/**********************************************************************************************************                      
     FACTORING ALGORTIHM using ECM (Elliptic Curve Method)                          
     =====================================================                      


     INPUT : Name of Parameterfile containg                                                                           

                 NOC :  Number Of Curves                                                                                         
                             v   :  >= second largest p | n                                                                          
                             w   :  Time spent on each curve                                                                     
                             n   :  The number to factorize                                                                        
                                                                                                                                                    
                                                                                                                                                    
     OUTPUT: List of numbers (Primes/Non Primes) Dividing n                                           
                                                                                                                                                   
     METHOD: Assume p | n, take point P=(p0 : p1 : 1) and hope                                             
                        P is a point on some elliptic curve. The calculate                                            
                        k _a P where k(v,w,n) is a great number. If we are                                       
                        lucky and p >= v, then we will find p or d : p | d | n.                                     
                        Then continue with n := d, n := n/d. For every n we,                                      
                        we try with NOC different a, P.                                                                        
                                                                                                                                                   
                                                                                                                                                   
                                                                                                                                                   
   __________________________                                                                                                     
   D. Jacquet and M. Skoglund                                                                                                  
******************************************************************************************************************/


const int ARRAY_SIZE=3;
const int PRIME_TEST_NUMBER=7;
const int BASE=10;
const int MAX_NUMBER_OF_DIGITS=1000;
const int MAX_NUMBER_OF_PRIMES=20000;
const int MAX_EXPONENT=1000000000;
const int NUMBER_OF_CURVES_IN_ONE_GROUP=10;  //-1 means all curves in one group

const int B2_START = 1000000000;
const int B2_STEP = 100000;


const float Probability_Reduction=0.10; 

struct List_Type{
  mpz_t val;
  struct List_Type *prev,*next;
};

// Read prime file
// =======================================================================
int calculate_number_of_primes(char *file_name){
  int number_of_primes,tmp;
  FILE *f;
  f=fopen(file_name,"r");
  if( f == NULL ){
    printf("Cannot open file %s \n",file_name);
    exit(99);
  }
  number_of_primes=0;
  while( !feof(f) ){
    tmp=0;
    fscanf(f,"%d",&tmp);
    if( tmp > 0 )
      number_of_primes++;
  }
  fclose(f);
  printf("The file %s has %d primes\n",file_name,number_of_primes);
  return number_of_primes;
}

int * read_prime_file(char *file_name, int number_of_primes){
  int i,*prime_array;
  FILE *f;
  f=fopen(file_name,"r");
  if( f == NULL ){
    printf("Cannot open file %s \n",file_name);
    exit(99);
  }
  prime_array = (int *)calloc(number_of_primes,sizeof(int));
  i=0;
  while( i<number_of_primes ){
    fscanf(f,"%d",&prime_array[i]);
    i++;
  }
  fclose(f);
  return prime_array;
}
//=======================================================================


// Managing List_Types
// =======================================================================
void Add_To_List(struct List_Type *list, mpz_t val){
  struct List_Type *tmp;
  
  if( list->prev == list && list->next == list && list != NULL){
    mpz_set(list->val,val);
    list->prev=list->next=NULL;
  }
  else if( list != NULL){
    tmp=list;
    while( tmp->next != NULL )
      tmp=tmp->next;
   
    tmp->next=(struct List_Type *)calloc(1,sizeof(struct List_Type));

    tmp->next->prev=tmp;
    mpz_set(tmp->next->val,val);
    tmp->next->next=NULL;
  }
  else{
    
    printf("Cannot Add to NULL -- list\n");
    exit(99);

  }
}
void Free_List(struct List_Type *list){
  struct List_Type *tmp,*next;
  tmp=list;
  while( tmp->prev != NULL )
    tmp=tmp->prev;
 
  while( tmp != NULL ){
    next=tmp->next;

    mpz_clear(tmp->val);
    free(tmp);
    
    tmp=next;    
  }
}
void Print_Factor_List(struct List_Type *Factor_List){
  struct List_Type *tmp,*start;

  start=Factor_List;
  while( start->prev != NULL )
    start=start->prev;
  printf("\n\n");
  printf("Algorithm Output\n");
  printf("================\n\n");
  printf("Prime Factors : ");
  
  tmp=start;
  while( tmp != NULL ){
    if( mpz_probab_prime_p(tmp->val,PRIME_TEST_NUMBER) >= 1 ){
      mpz_out_str(NULL,BASE,tmp->val);
      printf("\t");
    }
    tmp=tmp->next;
  }
  printf("\n");

  printf("Other Factors : ");
  tmp=start;
  while( tmp != NULL ){
    if( mpz_probab_prime_p(tmp->val,PRIME_TEST_NUMBER) == 0 ){
      mpz_out_str(NULL,BASE,tmp->val);
      printf("\t");
    }
    tmp=tmp->next;
  }
  printf("\n\n");

}
// =======================================================================


// mpz functions
//=======================================================================
int My_mpz_Log(mpz_t n, int Log_BASE){
  int exponent;
  mpz_t tmp;

  mpz_init(tmp);
  mpz_set_ui(tmp,Log_BASE);
  exponent=0;
  while( mpz_cmp(tmp,n) < 0 && exponent < MAX_EXPONENT ){
    
    mpz_mul_ui(tmp,tmp,Log_BASE);
    exponent++;
    
  }
  if( exponent == MAX_EXPONENT ){
    printf("You are trying to take log of a too big number, rewrite code...\n");
    exit(99);
  }
  mpz_clear(tmp);
  return exponent;
}
//=======================================================================



// Divide Out Start Primes
// =======================================================================
void Divide_Out(mpz_t n, mpz_t p,struct List_Type *Factor_List){
  mpz_t r;
  
  mpz_init(r);

  mpz_mod(r,n,p);
  while( mpz_cmp_ui(r,0) == 0 ){
    mpz_divexact(n,n,p);
    mpz_mod(r,n,p);
    
    Add_To_List(Factor_List,p);
  }

  mpz_clear(r);
}
void Divide_Out_With_Array(mpz_t n, int Number_Of_Primes, int *Prime_Array,struct List_Type *Factor_List){
  int i;
  mpz_t p;

  mpz_init(p);
  for( i=0 ; i<Number_Of_Primes ; i++){
    mpz_set_ui(p,Prime_Array[i]);
    Divide_Out(n,p,Factor_List);
  }
  mpz_clear(p);
}
// =======================================================================




// Array manipulating function
// =======================================================================
int Check_Zero(mpz_t *P){
  int i,ok=1;  

  if( mpz_cmp_ui(P[0],0) != 0 )
    ok=0;
  else if( mpz_cmp_ui(P[1],1) != 0 )
    ok=0;
  else if( mpz_cmp_ui(P[2],0) != 0 )
    ok=0;
  return ok;
}

void Set_Array(mpz_t *Target, mpz_t *Origin, int Size){
  int i;
  for( i=0 ; i<Size ; i++ )
    mpz_set(Target[i],Origin[i]);
}

void Free_Array(mpz_t *Array, int Size){
  int i;
  for( i=0 ; i<Size ; i++ )
    mpz_clear(Array[i]);
  free(Array);
}
void Free_Matrix(mpz_t **Matrix, int rows, int kols){
  int r,k;
  for( r=0 ; r<rows ; r++ ){
    for(k=0 ; k<kols ; k++)
      mpz_clear(Matrix[r][k]);
    free(Matrix[r]);
  }  
}
mpz_t * Get_Array(int Size){
  int i;
  mpz_t *Array;
  Array=(mpz_t *)calloc(Size,sizeof(mpz_t));
  for( i=0 ; i<Size ; i++ )
    mpz_init(Array[i]);
  return Array;
}

void Mod_Array(mpz_t *Array, mpz_t n, int Size){
  int i; 
  for( i=0 ; i<Size ; i++ )
    mpz_mod(Array[i],Array[i],n);
}
// =======================================================================


// ECM -- functions
// =======================================================================

// Get Random Values
// -----------------------------------------------------------------------
void Set_mpz_Random(mpz_t x,mpz_t n){
  int i,Number_Of_Digits=0;
  mpz_t tmp,tmp_BASE;

  mpz_init(tmp);
  mpz_init(tmp_BASE);
  mpz_set_ui(tmp,BASE);

  while( mpz_cmp(tmp,n) < 0 && Number_Of_Digits < MAX_NUMBER_OF_DIGITS ){
    
    mpz_mul_ui(tmp,tmp,BASE);
    Number_Of_Digits++;

  }

  mpz_set_ui(x,0);
  mpz_set_ui(tmp_BASE,1);

  for(i=0 ; i<Number_Of_Digits ; i++ ){

    mpz_mul_ui(tmp,tmp_BASE,rand() % BASE );
    mpz_mul_ui(tmp_BASE,tmp_BASE,BASE);

    mpz_add(x,x,tmp);

  }
  

  mpz_clear(tmp);
  mpz_clear(tmp_BASE);
}
void Set_Random_Value(mpz_t random, mpz_t n){
  mpz_set_ui(random,rand());
  mpz_mod(random,random,n);
}
mpz_t *Get_Random_Array(mpz_t n,int Size){
  int i;
  mpz_t *Array;
  Array=Get_Array(Size);
  for( i=0 ; i<Size ; i++)
    Set_Random_Value(Array[i],n);
  return Array;
}
mpz_t **Get_Random_Start_Points(mpz_t n,int Number_Of_Points){
  int i,j;
  mpz_t **Point_Array;

  Point_Array=(mpz_t **)calloc(Number_Of_Points,sizeof(mpz_t *));
  for(i=0 ; i<Number_Of_Points ; i++){
    Point_Array[i]=Get_Array(ARRAY_SIZE);
    for(j=0 ; j<ARRAY_SIZE-1 ; j++)
      Set_Random_Value(Point_Array[i][j],n);
    mpz_set_ui(Point_Array[i][ARRAY_SIZE-1],1);
  }
  return Point_Array;
}
// -----------------------------------------------------------------------


// Addition / Multiplication
// -----------------------------------------------------------------------
void Add(mpz_t d, mpz_t *R, mpz_t *P, mpz_t *Q,mpz_t a, mpz_t n){  
  int check; 
  mpz_t inv,lambda,tmp,r0,r1;
 
  mpz_set_ui(d,1);
  
  if( Check_Zero(P) == 1 )
    Set_Array(R,Q,ARRAY_SIZE);
  else if( Check_Zero(Q) == 1 )
    Set_Array(R,P,ARRAY_SIZE);
  else{
    
    mpz_init(inv);mpz_init(lambda);mpz_init(tmp);mpz_init(r0);mpz_init(r1);

    mpz_sub(tmp,P[0],Q[0]);
    check=mpz_invert(inv,tmp,n);
    
    if( check != 0 ){
 
      mpz_sub(lambda,P[1],Q[1]);
      mpz_mul(lambda,lambda,inv);
   
    } 
    else{
    
      mpz_gcd(d,tmp,n);
     
    }
    if( check == 0 && mpz_cmp(d,n) == 0 ){

      mpz_set_ui(d,1);
      mpz_add(tmp,P[1],Q[1]);
      check=mpz_invert(inv,tmp,n);
      
      if( check != 0 ){

	mpz_mul(lambda,P[0],P[0]);
	mpz_mul_ui(lambda,lambda,3);
	mpz_add(lambda,lambda,a);
	mpz_mul(lambda,lambda,inv);

      }
      else{
	mpz_gcd(d,tmp,n);	
      }
    }
    if( check != 0 ){

      mpz_mod(lambda,lambda,n);
      
      mpz_mul(r0,lambda,lambda);
      mpz_sub(r0,r0,P[0]);
      mpz_sub(r0,r0,Q[0]);
      
      mpz_sub(r1,P[0],r0);
      mpz_mul(r1,r1,lambda);
      mpz_sub(r1,r1,P[1]);
   
      mpz_set(R[0],r0);
      mpz_set(R[1],r1);
    }
    
    mpz_clear(inv);mpz_clear(lambda);mpz_clear(tmp);mpz_clear(r0);mpz_clear(r1);
    
  }
  
  Mod_Array(R,n,ARRAY_SIZE);

}
void Mul_Under(mpz_t d,mpz_t k, mpz_t *R, mpz_t *P,mpz_t a, mpz_t n){
  mpz_t tmp_k,tmp_2,*tmp_P;

  tmp_P=(mpz_t *)calloc(ARRAY_SIZE,sizeof(mpz_t));
  Set_Array(tmp_P,P,ARRAY_SIZE);
  mpz_init(tmp_k);mpz_init(tmp_2);

  mpz_set_ui(tmp_k,2);
  mpz_set_ui(tmp_2,2);
  mpz_set_ui(d,1);

  while( mpz_cmp(tmp_k,k) < 0 && mpz_cmp_ui(d,1) == 0 ){
    Add(d,tmp_P,tmp_P,tmp_P,a,n);
    mpz_mul_ui(tmp_k,tmp_k,2);
  }
  mpz_divexact(tmp_k,tmp_k,tmp_2);
  if( mpz_cmp_ui(d,1) == 0 )
    Add(d,R,R,tmp_P,a,n);
  Free_Array(tmp_P,ARRAY_SIZE);
  mpz_sub(k,k,tmp_k);
  mpz_clear(tmp_k);
  mpz_clear(tmp_2);
}
void Mul(mpz_t d,mpz_t *R, mpz_t k, mpz_t *P,mpz_t a, mpz_t n){

  mpz_t tmp_k,k_diff,index;

  mpz_init(tmp_k);
  mpz_init(index);

  mpz_set(tmp_k,k);
  mpz_set_ui(d,1);
  mpz_set_ui(index,0);

  mpz_sub_ui(tmp_k,tmp_k,1);
  Set_Array(R,P,ARRAY_SIZE);

  while( mpz_cmp_ui(tmp_k,0) > 0 && mpz_cmp_ui(d,1) == 0 ){

    Mul_Under(d,tmp_k,R,P,a,n);
    mpz_add_ui(index,index,1);

  }

  if( mpz_cmp_ui(d,1) > 0 && mpz_cmp(d,n) ){
    mpz_init(k_diff);
    mpz_sub(k_diff,k,tmp_k);

    mpz_out_str(NULL,BASE,k_diff);printf(" , ");
    mpz_out_str(NULL,BASE,index);printf(" ");

    mpz_clear(k_diff);
  }

  mpz_clear(tmp_k);
  mpz_clear(index);
}
// -----------------------------------------------------------------------


// Calculate parameters
//-----------------------------------------------------------------------
int *Calculate_E_Array(int B1_index, int *Prime_Array, mpz_t C){
  int r;
  int *E_Array;
  E_Array = (int *)calloc(B1_index,sizeof(int));
  for( r=0 ; r<B1_index ; r++){
    E_Array[r]=My_mpz_Log(C,Prime_Array[r]);
  }
  return E_Array;
}
float Probability_With_One_Curve(mpz_t C, int B1){
  float prob;
  int tmp_Log_BASE;
  int tmp_int_Exponent;
  float tmp_float_Exponent;
  
  tmp_Log_BASE=2;
  tmp_int_Exponent=My_mpz_Log(C,tmp_Log_BASE);
  tmp_float_Exponent=log(log(tmp_Log_BASE+0.0)*(tmp_int_Exponent+0.0));
  tmp_float_Exponent = (1.0-tmp_float_Exponent) / log(B1);

  prob = exp(tmp_float_Exponent * (log(tmp_Log_BASE+0.0)*(My_mpz_Log(C,tmp_Log_BASE) + 0.0)) );
  
  return prob;
}
void Calculate_C(mpz_t C, mpz_t B2){
  mpz_sqrt(C,B2);
  mpz_sub(C,B2,C);
  mpz_add_ui(C,C,1);
}
int Calculate_B1_index(mpz_t C,int *E_Array, int *Prime_Array, int Size){
  int r,index,Log_Base;
  float work,value, max_value,sum,total_sum;
  Log_Base=2;
  max_value=-1;
  work=0.0;
  index=0;
  sum=total_sum=0;

  for(r=0 ; r<Size ; r++)
    total_sum += log(Prime_Array[r])*(E_Array[r]+0.0)*log(Log_Base+0.0);

  for(r=0 ; r<Size ; r++){
    sum += log(Prime_Array[r])*(E_Array[r]+0.0)*log(Log_Base+0.0);
    value=total_sum*Probability_With_One_Curve(C,Prime_Array[r])/sum;
    if( value > max_value ){
      index=r;
      max_value=value;
    }
  }
  return index;
}
int Calculate_NOC(mpz_t C, int B1 ){
  float NOC_float,prob;
  int NOC;
  prob=Probability_With_One_Curve(C,B1);
  printf("prob=%f\n",prob);
  NOC_float=log(1.0-Probability_Reduction)/log(1.0-prob);
  NOC=NOC_float;
  return NOC;
}
//-----------------------------------------------------------------------


// Factor Function
//-----------------------------------------------------------------------
void Factor_With_Parameters(mpz_t d, mpz_t n, int B1_index, mpz_t B2, int NOC, int *E_Array, struct List_Type *Factor_List){
  int j,r,curve,tmp_curve,curve_group;
  int Number_Of_Groups;
  int divided;
  mpz_t mpz_r;
  mpz_t **P;
  mpz_t *a;
  
  mpz_init(mpz_r);

  if( NUMBER_OF_CURVES_IN_ONE_GROUP == -1 )
    Number_Of_Groups=1;
  else
    Number_Of_Groups= (NOC-1) / NUMBER_OF_CURVES_IN_ONE_GROUP + 1;

  a=Get_Random_Array(n,NOC);
  P=Get_Random_Start_Points(n,NOC);
  divided=0;
  curve_group=0;
  while( curve_group<Number_Of_Groups && divided == 0){
    curve= NUMBER_OF_CURVES_IN_ONE_GROUP*curve_group;
    r=0;
    while( r<B1_index && divided == 0){
      tmp_curve=0;
      while( curve+tmp_curve < NOC && (Number_Of_Groups == 1  || tmp_curve<NUMBER_OF_CURVES_IN_ONE_GROUP)  && divided == 0){
	j=0;
	while( j < E_Array[r] ){ 
	
	  mpz_set_ui(mpz_r,r);
	  Mul(d,P[curve+tmp_curve],mpz_r,P[curve+tmp_curve],a[curve+tmp_curve], n);
	  if( mpz_cmp_ui(d,1) > 0 && mpz_cmp(d,n) != 0  ){
	    divided=1;
	  }
	j++;
	}
      tmp_curve++;
      }
      r++;
    }
    curve_group++;
  }
  mpz_clear(mpz_r);
  Free_Array(a,NOC);
  Free_Matrix(P,NOC,ARRAY_SIZE);
}
void Factor_Integer(mpz_t n,mpz_t B2,int Size_Of_Prime_Array,int *Prime_Array,struct List_Type *Factor_List){
  int B1_index,NOC,*E_Array,Size;
  mpz_t d,new_n,C,sqrt_n;

  if ( mpz_probab_prime_p(n,PRIME_TEST_NUMBER) == 0  && mpz_cmp_ui(n,1) > 0 ){

    mpz_init(d);
    mpz_init(new_n);
    mpz_init(sqrt_n);
    mpz_init(C); 
    
    mpz_sqrt(sqrt_n,n);
    
    Calculate_C(C,B2);
    E_Array=Calculate_E_Array(Size_Of_Prime_Array,Prime_Array,C);
    B1_index=Calculate_B1_index(C,E_Array,Prime_Array,Size_Of_Prime_Array);
    NOC=Calculate_NOC(C, Prime_Array[B1_index] );
    
    printf("Factoring ");mpz_out_str(NULL,BASE,n);printf(" with ");printf(" B1_index %7d B1 %12d NOC %12d",B1_index,Prime_Array[B1_index],NOC);printf(" B2 ");mpz_out_str(NULL,BASE,B2);printf("\n");
    
    Factor_With_Parameters(d,n,B1_index,B2,NOC,E_Array,Factor_List);
    if( mpz_cmp_ui(d,1) > 0 ){
      mpz_divexact(new_n,n,d);
      mpz_set_ui(B2,B2_START);
      Factor_Integer(d,B2,Size_Of_Prime_Array,Prime_Array,Factor_List);
      Factor_Integer(new_n,B2,Size_Of_Prime_Array,Prime_Array,Factor_List);
    }
    else{
      if( mpz_cmp(B2,sqrt_n) <= 0 ){
	mpz_mul_ui(B2,B2,B2_STEP);
      }
      Factor_Integer(n,B2,Size_Of_Prime_Array,Prime_Array,Factor_List);
    }
    free(E_Array);
    mpz_clear(d);
    mpz_clear(B2);
    mpz_clear(new_n);
    mpz_clear(sqrt_n);
    mpz_clear(C);
  }
  else{
    if(mpz_cmp_ui(n,1) > 0)
      Add_To_List(Factor_List,n);
  }
}
//-----------------------------------------------------------------------



int main(int Number_Of_Arguments, char **Argument_Array){
  mpz_t n,B2;
  int number_of_primes;
  int *prime_array;
  char prime_file[80]="/home/matte/martins/matte/primtal.lst";
  struct List_Type *Factor_List;

  if( Number_Of_Arguments != 2 ){
    printf("You must type in one (1) argument, the number to factor. You have entered %d argumnets\n",Number_Of_Arguments-1);
    exit(99);
  }
  srandom(time(0));

  Factor_List=(struct List_Type *)calloc(1,sizeof(struct List_Type));
  Factor_List->prev=Factor_List->next=Factor_List;
  
  number_of_primes=calculate_number_of_primes(prime_file);
  if(number_of_primes > MAX_NUMBER_OF_PRIMES )
    number_of_primes=MAX_NUMBER_OF_PRIMES;  
  prime_array=read_prime_file(prime_file,number_of_primes);

  mpz_init(B2);
  mpz_set_ui(B2,B2_START);
  
  mpz_init(n);
  mpz_set_str(n,Argument_Array[1],BASE);

  Divide_Out_With_Array(n,number_of_primes,prime_array,Factor_List);
  Factor_Integer(n,B2,number_of_primes,prime_array,Factor_List);

  Print_Factor_List(Factor_List);
  Free_List(Factor_List);

  mpz_clear(B2);
  mpz_clear(n);
  return 0;
}
