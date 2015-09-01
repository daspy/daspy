'''
Copyright of DasPy:
Author - Xujun Han (Forschungszentrum Jülich, Germany)
x.han@fz-juelich.de, xujunhan@gmail.com

DasPy was funded by:
1. Forschungszentrum Jülich, Agrosphere (IBG 3), Jülich, Germany
2. Cold and Arid Regions Environmental and Engineering Research Institute, Chinese Academy of Sciences, Lanzhou, PR China
3. Centre for High-Performance Scientific Computing in Terrestrial Systems: HPSC TerrSys, Geoverbund ABC/J, Jülich, Germany

Please include the following references related to DasPy:
1. Han, X., Li, X., He, G., Kumbhar, P., Montzka, C., Kollet, S., Miyoshi, T., Rosolem, R., Zhang, Y., Vereecken, H., and Franssen, H. J. H.: 
DasPy 1.0 : the Open Source Multivariate Land Data Assimilation Framework in combination with the Community Land Model 4.5, Geosci. Model Dev. Discuss., 8, 7395-7444, 2015.
2. Han, X., Franssen, H. J. H., Rosolem, R., Jin, R., Li, X., and Vereecken, H.: 
Correction of systematic model forcing bias of CLM using assimilation of cosmic-ray Neutrons and land surface temperature: a study in the Heihe Catchment, China, Hydrology and Earth System Sciences, 19, 615-629, 2015a.
3. Han, X., Franssen, H. J. H., Montzka, C., and Vereecken, H.: 
Soil moisture and soil properties estimation in the Community Land Model with synthetic brightness temperature observations, Water Resour Res, 50, 6081-6105, 2014a.
4. Han, X., Franssen, H. J. H., Li, X., Zhang, Y. L., Montzka, C., and Vereecken, H.: 
Joint Assimilation of Surface Temperature and L-Band Microwave Brightness Temperature in Land Data Assimilation, Vadose Zone J, 12, 0, 2013.
'''

import os,socket,numpy,scipy.weave


def ParFor_Common(CLM_NA,ColumnNum,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,Array_Origin,varin,Mask_Region,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue):
    
    varout = numpy.zeros(ColumnNum,dtype=numpy.float32)
    
    support_code="""
        #define BZ_THREADSAFE
        #include <blitz/array.h>
    """
    code="""
    
        using namespace std;
        using namespace blitz;
        
        int Column_Index,  Row_Index, Col_Index;
        
        omp_set_num_threads(omp_get_num_procs_ParFor);
        
        #pragma omp parallel private(Column_Index,Row_Index,Col_Index) shared(varout, varin)
        {
            #pragma omp for schedule(static)
            for(Column_Index=0;Column_Index<ColumnNum;Column_Index++){
                cout<<flush;
                Row_Index=Row_Numbers-cols1d_jxy(Column_Index);
                Col_Index=cols1d_ixy(Column_Index) - 1;
                varout(Column_Index)=Array_Origin(Column_Index);
                //cout<<varin(Row_Index,Col_Index)<<" "<<Row_Index<<" "<<Col_Index<<endl;
                if(!Mask_Region(Row_Index,Col_Index)){
                    varout(Column_Index)=varin(Row_Index,Col_Index);
                }
            }
        }
    """
    
    args =['-pthread -O3 -fopenmp']
    if socket.gethostname()[0:4] == 'node':
        args =['-pthread -O3 -fopenmp']
    if os.name == 'nt' or os.name == 'posix':
        libs=['gomp']
    if os.name == 'posix' and socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if os.name == 'posix' and socket.gethostname()[0:4] == 'node':
        libs=['gomp','blas','lapack']
    dirs1=['.',DAS_Depends_Path + 'include',DAS_Depends_Path + 'include/python2.7']
    dirs2=['.',DAS_Depends_Path + 'lib',DAS_Depends_Path + 'lib64']
    vars_list = ['CLM_NA','ColumnNum','Row_Numbers','Col_Numbers','cols1d_ixy','cols1d_jxy','Array_Origin','varin','Mask_Region','varout','omp_get_num_procs_ParFor','NAvalue']
    scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])
    if socket.gethostname()[0:4] == 'node':
        scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])

    #print numpy.min(varout),numpy.max(varout)
    return varout

def ParFor_Fusion(Def_ReBEL,Obs_Grid_Length,CLM_NA,Ensemble_Number,Row_Numbers,Col_Numbers,Obs_Index_Dim,Mask_Index,Def_Localization,num_local_obs, Mask,pos_obs,loc_len,loc_function,eps,kappa,Obs_Grid,Observation_Matrix,\
                  Prop_Grid_Array_Sys_Ens_Mean,Prop_Grid_Array_H_Trans_Ens_Mean,Prop_Grid_Array_Sys,Prop_Grid_Array_H_Trans,K,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue,Def_Print):
    
    Analysis_Grid = numpy.zeros((Row_Numbers,Col_Numbers),dtype=numpy.float32)
    Analysis_Grid_Array = numpy.zeros((Ensemble_Number,Row_Numbers,Col_Numbers),dtype=numpy.float32)
    
    support_code="""
        #include <armadillo>
        #include <cmath>
        #include <cstring>
        #include <boost/math/special_functions/bessel.hpp>
        #include <boost/math/special_functions/gamma.hpp>
        using namespace arma;
        using namespace std;
        using namespace boost;
        
        //#ifndef isnan
        //inline bool isnan(double x) {
        //    return x != x;
        //}
        //#endif
        
        vec Exponential(vec h, double r) {
            int n = h.n_elem;
            vec corelation(n);
            for (int i = 0; i < n; i++) {
                if (h(i) == 0.0) {
                    corelation(i) = 1.0;
                } else {
                    corelation(i) = exp(-h(i) / r);
                }
            }
            return corelation;
        }
        
        vec Gaussian(vec h, double r) {
            int n = h.n_elem;
            vec corelation(n);
            for (int i = 0; i < n; i++) {
                if (h(i) == 0.0) {
                    corelation(i) = 1.0;
                } else {
                    double hr = h(i) / r;
                    corelation(i) = exp(-(hr * hr));
                }
            }
            return corelation;
        }
        
        vec Spherical(vec h, double r) {
            int n = h.n_elem;
            vec corelation(n);
            for (int i = 0; i < n; i++) {
                if (h(i) == 0.0) {
                    corelation(i) = 1.0;
                }
                if (h(i) >= r) {
                    corelation(i) = 0.0;
                } else {
                    double hr = h(i) / r;
                    corelation(i) = 1.0 - hr * (1.5 - 0.5 * hr * hr);
                }
            }
            return corelation;
        }
        
        vec SteMat(vec h, double r, double kappa) {
            int n = h.n_elem;
            vec corelation(n);
            for (int i = 0; i < n; i++) {
                if (h(i) == 0.0) {
                    corelation(i) = 1.0;
        
                } else {
                    double hr = h(i) / r;
                    double bes = boost::math::cyl_bessel_k(kappa, 2.0 * sqrt(kappa) * hr);
                    if (!isfinite(bes)) {
                        corelation(i) = 1.0;
                    }
                    if (bes == 0.0) {
                        corelation(i) = 0.0;
                    } else {
                        double mult = pow(2.0, (1.0 - kappa)) / boost::math::tgamma(kappa) * pow(
                                (2.0 * sqrt(kappa) * hr), kappa);
                        if (!isfinite(mult)) {
                            corelation(i) = 0.0;
                        } else {
                            corelation(i) = bes * mult;
                        }
                    }
                }
        
            }
            return corelation;
        }
        
        vec Matern(vec h, double r, double kappa) {
            int n = h.n_elem;
            vec corelation(n);
            for (int i = 0; i < n; i++) {
                if (h(i) == 0.0) {
                    corelation(i) = 1.0;
                }
                if (h(i) > 600 * r) {
                    corelation(i) = 0.0;
                } else {
                    double hr = h(i) / r;
                    double r1 = pow(2, (kappa - 1.0)) * boost::math::tgamma(kappa);
                    double bes = boost::math::cyl_bessel_k(kappa, hr);
                    corelation(i) = (1.0 / r1) * pow(hr, kappa) * bes;
                }
            }
            return corelation;
        }
        
        // r is the range or distance parameter (r>0) which measures how quickly the correlations decay with distance
        // v is the smoothness parameter (v>0)
        // if v moves towards to infinity, the Matern function corresponds to a Gaussian model
        // if v is equal to 0.5, the Matern function is Exponential
    
        vec calc_loccoeffs(double radius, string tag, vec dist, double kappa)
        {
            vec coeffs(dist.n_elem);
            double R;
            
            if(tag == "Gauss"){
                R = radius;
                coeffs = exp(-0.5 * pow((dist / R) , 2));
                //cout<<R<<" "<<coeffs<<endl;
            }
            else if( tag == "Gaspari_Cohn"){
                R = radius * 1.7386;
                uvec ind1 = find(dist <= R);
                vec r2 = pow((dist.elem(ind1) / R), 2);
                vec r3 = pow((dist.elem(ind1) / R), 3);
                coeffs.elem(ind1) = 1 + r2 * (- r3 / 4.0 + r2 / 2.0) + r3 * (5.0 / 8.0) - r2 * (5.0 / 3.0);
                uvec ind2_temp1 = find(dist > R);
                uvec ind2_temp2 = find(dist <= (R * 2.0));
                uvec ind2 =  ind2_temp1.elem(find(ind2_temp1 == ind2_temp2));
                vec r1 = (dist.elem(ind2) / R);
                r2 = pow((dist.elem(ind2) / R) , 2);
                r3 = pow((dist.elem(ind2) / R) , 3);
                coeffs.elem(ind2) = r2 * (r3 / 12.0 - r2 / 2.0) + r3 * (5.0 / 8.0) + r2 * (5.0 / 3.0) - r1 * 5.0 + 4.0 - (2.0 / 3.0) / r1;
            }
            else if( tag == "Cosine"){
                R = radius * 2.3167;
                //ind = numpy.nonzero(dist <= R)[0];
                //r = dist[ind] / R;
                //coeffs[ind] = (1.0 + numpy.cos(r * numpy.pi)) / 2.0;
            }
            else if( tag == "Cosine_Squared"){
                R = radius * 3.2080;
                //ind = numpy.nonzero(dist <= R)[0];
                //r = dist[ind] / R;
                //coeffs[ind] = ((1 + numpy.cos(r * numpy.pi)) / 2.0) ** 2;
            }
            else if( tag == "Lewitt"){
            
                R = radius * 4.5330;
                //m = 10.0;
                //alpha = 1.0;
                //ind = numpy.nonzero(dist <= R)[0];
                //for i in ind:
                  //  r = dist(i) / R;
                    //rr = 1.0 - r ** 2.0;
                    //coeffs(i) = rr ** (m/2.0) * besseli(m, alpha * rr ** 0.5) / besseli(m, alpha);
            }
            else if( tag == "Exp3"){
        
                R = radius;
                //coeffs = numpy.exp(-0.5 * (dist / R) ** 3.0);
            }
            else if( tag == "Cubic"){
                
                R = radius * 1.8676;
                //ind = numpy.nonzero(dist < R)[0];
                //coeffs[ind] = (1.0 - (dist[ind] / R) ** 3.0) ** 3.0;
            }
            else if( tag == "Quadro"){
                
                R = radius * 1.7080;
                //ind = numpy.nonzero(dist < R)[0];
                //coeffs[ind] = (1.0 - (dist[ind] / R) ** 4.0) ** 4.0;
            }
            else if( tag == "Step"){
              
                R = radius;
                //ind = numpy.nonzero(dist < R)[0];
                //coeffs[ind] = 1.0;
            }
            else if( tag == "SteMat"){
                R = radius;
                coeffs = SteMat(dist,R,kappa);
            }
            else if( tag == "Matern"){
                  R = radius;
                  coeffs = Matern(dist,R,kappa);
            }
            else if( tag == "Exponential"){
                R = radius;
                coeffs = Exponential(dist,R);
            }
            else if( tag == "Gaussian"){
                R = radius;
                coeffs = Gaussian(dist,R);
            }
            else if( tag == "Spherical"){
                R = radius;
                coeffs = Spherical(dist,R);
            }
            else{
                cout<<"Wrong coorelation type name!!!"<<endl;
            }
            return coeffs;
        }
        
        void find_localobs(double loc_len, string loc_function, vec pos_sys, mat pos_obs, double eps, double kappa, uvec &obs_index, vec &coeffs)
        {
            // 2D case
            double ii = pos_sys(0);
            double jj = pos_sys(1);
            
            vec dist = sqrt(pow(abs(ii - pos_obs.col(0)), 2) + pow(abs(jj - pos_obs.col(1)), 2));
            //cout<<dist<<endl;
            vec coeffs_temp = calc_loccoeffs(loc_len, loc_function, dist, kappa);
            coeffs_temp = coeffs_temp / max(1.0,coeffs_temp.max());
            //cout<<coeffs_temp<<endl;
            //print dist
            //print coeffs
            
            uvec obs_index_temp = find(coeffs_temp >= eps);
            //cout<<obs_index_temp<<endl;
            vec coeffs_temp2 = coeffs_temp.elem(obs_index_temp);
            //cout<<coeffs_temp2<<endl;
            
            if(coeffs_temp2.n_elem > 4){
                uvec indices = sort_index(coeffs_temp2,1);
                coeffs = coeffs_temp2.elem(indices.subvec(0,3));
                obs_index = obs_index_temp.elem(indices.subvec(0,3));
                //cout<<coeffs<<" "<<obs_index<<endl;
            }
            else{
                coeffs = coeffs_temp2;
                obs_index = obs_index_temp;
            }
            //cout<<coeffs<<" "<<obs_index<<endl;
            //cout<<"******************************"<<endl;
            //cout<<coeffs<<endl;
        }
    """
    code="""        
        using namespace std;
        using namespace arma;
            
        int Obs_Index, Row_Index, Col_Index, Ens_Index, num_local_obs_temp;
        double Observation_Temp;
        uvec obs_index_vec;
        vec coeffs,pos_sys(2);
        mat Obs_Grid_Trans(Obs_Grid_Length,3), pos_obs_Trans(Obs_Grid_Length,2);
        vec Obs_Grid_Col;
        
        omp_set_num_threads(omp_get_num_procs_ParFor);
        
         #pragma omp parallel private(Row_Index,Col_Index) shared(Obs_Grid_Trans,pos_obs_Trans)
        {
            #pragma omp for schedule(static) collapse(2)
            for(Row_Index=0; Row_Index<Obs_Grid_Length;Row_Index++){
                for(Col_Index=0; Col_Index<3;Col_Index++){
                    //cout<<Row_Index<<" "<<Col_Index<<endl;
                    cout<<flush;
                    Obs_Grid_Trans(Row_Index, Col_Index) = Obs_Grid(Row_Index, Col_Index);
                    if(Col_Index <=1 ){
                        pos_obs_Trans(Row_Index, Col_Index) = pos_obs(Row_Index, Col_Index);
                    }
                }
            }
        }
        
        #pragma omp parallel private(Obs_Index,Row_Index,Col_Index,Ens_Index,pos_sys,coeffs,num_local_obs_temp,Observation_Temp,Obs_Grid_Col,obs_index_vec) shared(Analysis_Grid,Analysis_Grid_Array)
        {
            #pragma omp for schedule(static) collapse(2)
            
            for(Row_Index=0; Row_Index<Row_Numbers;Row_Index++){
                for(Col_Index=0; Col_Index<Col_Numbers;Col_Index++){
                    //cout<<Row_Index<<" "<<Col_Index<<endl;
                    cout<<flush;
                    pos_sys.set_size(2);
                    if(Mask_Index(Col_Index + Row_Index * Col_Numbers)){
                        Analysis_Grid(Row_Index, Col_Index) = NAvalue;  // Apply the Boundary Mask
                        for(Ens_Index=0;Ens_Index<Ensemble_Number;Ens_Index++){
                            Analysis_Grid_Array(Ens_Index, Row_Index, Col_Index) = NAvalue;
                        }
                    }
                    else{
                        if(Def_Localization){
                            pos_sys(0) = Mask(Col_Index + Row_Index * Col_Numbers, 0);
                            pos_sys(1) = Mask(Col_Index + Row_Index * Col_Numbers, 1);
                            //cout<<pos_sys(0)<<" "<<pos_sys(1)<<endl;
                            #pragma omp critical
                            {
                                find_localobs(loc_len, loc_function, pos_sys, pos_obs_Trans, eps, kappa, obs_index_vec, coeffs);
                            }
                            if (Def_Print)
                                cout<<"There are "<<obs_index_vec.n_elem<<" Observations"<<endl;
                            if(obs_index_vec.n_elem > 0){
                                //coeffs_max_index = numpy.argmax(coeffs);
                                //Observation_Temp = Obs_Grid_Trans(obs_index_vec(coeffs_max_index),2);
                                //Observation_Temp = mean(Obs_Grid_Trans.col(2).elem(obs_index_vec));
                                num_local_obs_temp = min(int(num_local_obs),int(obs_index_vec.n_elem));
                                coeffs = (coeffs / max(coeffs));
                                Obs_Grid_Col = Obs_Grid_Trans.col(2);
                                Observation_Temp = sum(sum(Obs_Grid_Col.elem(obs_index_vec.subvec(0,num_local_obs_temp-1))%coeffs.subvec(0,num_local_obs_temp-1)/sum(coeffs.subvec(0,num_local_obs_temp-1))));
                            }
                            else{
                                Observation_Temp = Observation_Matrix(Row_Index, Col_Index);
                            }
                            coeffs.zeros();
                            obs_index_vec.zeros();
                            //cout<<Observation_Temp<<endl;
                            //Observation_Temp = Observation_Matrix(Row_Index, Col_Index);
                        }
                        else{
                            Observation_Temp = Observation_Matrix(Row_Index, Col_Index);
                        }
                        if(Observation_Temp != NAvalue && Prop_Grid_Array_Sys_Ens_Mean(Row_Index, Col_Index) != NAvalue){
                            if (Def_ReBEL == 0){
                                Analysis_Grid(Row_Index, Col_Index) = Prop_Grid_Array_Sys_Ens_Mean(Row_Index, Col_Index) + K * (Observation_Temp - Prop_Grid_Array_H_Trans_Ens_Mean(Row_Index, Col_Index));
                                for(Ens_Index=0;Ens_Index<Ensemble_Number;Ens_Index++){
                                    Analysis_Grid_Array(Ens_Index, Row_Index, Col_Index) = Prop_Grid_Array_Sys(Ens_Index, Row_Index, Col_Index) + K * (Observation_Temp - Prop_Grid_Array_H_Trans(Ens_Index, Row_Index, Col_Index));
                                }
                            }
                            else if (Def_ReBEL == 2){
                                Analysis_Grid(Row_Index, Col_Index) = Observation_Temp;
                                for(Ens_Index=0;Ens_Index<Ensemble_Number;Ens_Index++){
                                    Analysis_Grid_Array(Ens_Index, Row_Index, Col_Index) = Observation_Temp;
                                }
                            }
                        }
                        else{
                            Analysis_Grid(Row_Index, Col_Index) = Prop_Grid_Array_Sys_Ens_Mean(Row_Index, Col_Index);
                            for(Ens_Index=0;Ens_Index<Ensemble_Number;Ens_Index++){
                                Analysis_Grid_Array(Ens_Index, Row_Index, Col_Index) = Prop_Grid_Array_Sys(Ens_Index, Row_Index, Col_Index);
                            }
                        }
                        
                        //cout<<"Analysis is "<<Analysis_Grid(Row_Index, Col_Index)<<endl;
                        if(isnan(Analysis_Grid(Row_Index,Col_Index))){
                            cout<< Observation_Temp <<endl;
                        }
                    }
                }
            }
        }
    """
    
    args =['-pthread -O3 -fopenmp']
    if socket.gethostname()[0:4] == 'node':
        args =['-pthread -O3 -fopenmp']
    if os.name == 'nt' or os.name == 'posix':
        libs=['gomp','armadillo']
    if os.name == 'posix' and socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack','armadillo']
    if socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack','armadillo']
    if os.name == 'posix' and socket.gethostname()[0:4] == 'node':
        libs=['gomp','blas','lapack','armadillo']
    dirs1=['.',DAS_Depends_Path + 'include',DAS_Depends_Path + 'include/python2.7']
    dirs2=['.',DAS_Depends_Path + 'lib',DAS_Depends_Path + 'lib64']
    vars_list = ['Def_ReBEL','Obs_Grid_Length','CLM_NA','Ensemble_Number','Row_Numbers','Col_Numbers','Obs_Index_Dim','Mask_Index','Def_Localization','num_local_obs','Mask','pos_obs','loc_len','loc_function','eps','kappa','Obs_Grid','Observation_Matrix',\
            'Prop_Grid_Array_Sys_Ens_Mean','Prop_Grid_Array_H_Trans_Ens_Mean','Prop_Grid_Array_Sys','Prop_Grid_Array_H_Trans','K','Analysis_Grid','Analysis_Grid_Array','omp_get_num_procs_ParFor','NAvalue','Def_Print']
    scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])
    if socket.gethostname()[0:4] == 'node':
        scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])

    return Analysis_Grid, Analysis_Grid_Array


def ParFor_PFT_Prepare_Model(Row_Numbers,Col_Numbers,maxpft,Parameter_PFT_Space_Ensemble,PFT_Dominant_Index,Analysis_Grid_Matrix,DAS_Depends_Path,omp_get_num_procs_ParFor):
        
    support_code="""
        #define BZ_THREADSAFE
        #include <blitz/array.h>
    """
    code="""
    
        using namespace std;
        using namespace blitz;
        
        int PFT_Type_Index, Row_Index, Col_Index;
        
        omp_set_num_threads(omp_get_num_procs_ParFor);
        
        #pragma omp parallel private(PFT_Type_Index,Row_Index,Col_Index) shared(Analysis_Grid_Matrix,Parameter_PFT_Space_Ensemble)
        {
            #pragma omp for schedule(static) collapse(2)
            for(Row_Index=0;Row_Index<Row_Numbers;Row_Index++){
                for(Col_Index=0;Col_Index<Col_Numbers;Col_Index++){
                    PFT_Type_Index = int(PFT_Dominant_Index(Row_Index,Col_Index));
                    Analysis_Grid_Matrix(Row_Index,Col_Index) = Parameter_PFT_Space_Ensemble(PFT_Type_Index,Row_Index,Col_Index);
                }
            }
        }
    """
    
    args =['-pthread -O3 -fopenmp']
    if socket.gethostname()[0:4] == 'node':
        args =['-pthread -O3 -fopenmp']
    if os.name == 'nt' or os.name == 'posix':
        libs=['gomp']
    if os.name == 'posix' and socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if os.name == 'posix' and socket.gethostname()[0:4] == 'node':
        libs=['gomp','blas','lapack']
    dirs1=['.',DAS_Depends_Path + 'include',DAS_Depends_Path + 'include/python2.7']
    dirs2=['.',DAS_Depends_Path + 'lib',DAS_Depends_Path + 'lib64']
    vars_list = ['Row_Numbers','Col_Numbers','maxpft','Parameter_PFT_Space_Ensemble','PFT_Dominant_Index','Analysis_Grid_Matrix','omp_get_num_procs_ParFor']
    scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])
    if socket.gethostname()[0:4] == 'node':
        scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])

    #print numpy.min(Parameter_PFT_Space_Ensemble),numpy.max(Parameter_PFT_Space_Ensemble)
    return

def ParFor_PFT(Row_Numbers,Col_Numbers,maxpft,Parameter_PFT_Space_Ensemble,PFT_Dominant_Index,Analysis_Grid_Matrix,DAS_Depends_Path,omp_get_num_procs_ParFor):
        
    support_code="""
        #define BZ_THREADSAFE
        #include <blitz/array.h>
    """
    code="""
    
        using namespace std;
        using namespace blitz;
        
        int PFT_Type_Index, Row_Index, Col_Index;
        
        omp_set_num_threads(omp_get_num_procs_ParFor);
        
        #pragma omp parallel private(PFT_Type_Index,Row_Index,Col_Index) shared(Parameter_PFT_Space_Ensemble)
        {
            #pragma omp for schedule(static) collapse(2)
            for(Row_Index=0;Row_Index<Row_Numbers;Row_Index++){
                for(Col_Index=0;Col_Index<Col_Numbers;Col_Index++){
                    PFT_Type_Index = int(PFT_Dominant_Index(Row_Index,Col_Index));
                    Parameter_PFT_Space_Ensemble(PFT_Type_Index,Row_Index,Col_Index) = Analysis_Grid_Matrix(Row_Index,Col_Index);
                }
            }
        }
    """
    
    args =['-pthread -O3 -fopenmp']
    if socket.gethostname()[0:4] == 'node':
        args =['-pthread -O3 -fopenmp']
    if os.name == 'nt' or os.name == 'posix':
        libs=['gomp']
    if os.name == 'posix' and socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if os.name == 'posix' and socket.gethostname()[0:4] == 'node':
        libs=['gomp','blas','lapack']
    dirs1=['.',DAS_Depends_Path + 'include',DAS_Depends_Path + 'include/python2.7']
    dirs2=['.',DAS_Depends_Path + 'lib',DAS_Depends_Path + 'lib64']
    vars_list = ['Row_Numbers','Col_Numbers','maxpft','Parameter_PFT_Space_Ensemble','PFT_Dominant_Index','Analysis_Grid_Matrix','omp_get_num_procs_ParFor']
    scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])
    if socket.gethostname()[0:4] == 'node':
        scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])

    #print numpy.min(Parameter_PFT_Space_Ensemble),numpy.max(Parameter_PFT_Space_Ensemble)
    return


    
def ParFor_PFT_Block_Assim(Row_Numbers,Col_Numbers,maxpft,Parameter_Space_SubBlock,Parameter_PFT_Space_Ensemble_SubBlock_Sens,PFT_Dominant_Index,DAS_Depends_Path,omp_get_num_procs_ParFor):
        
    support_code="""
        #define BZ_THREADSAFE
        #include <blitz/array.h>
    """
    code="""
    
        using namespace std;
        using namespace blitz;
        
        int PFT_Type_Index, Row_Index, Col_Index;
        
        omp_set_num_threads(omp_get_num_procs_ParFor);
        
        #pragma omp parallel for schedule(static) collapse(2) private(PFT_Type_Index,Row_Index,Col_Index) shared(Parameter_Space_SubBlock)
        for(Row_Index=0;Row_Index<Row_Numbers;Row_Index++){
            for(Col_Index=0;Col_Index<Col_Numbers;Col_Index++){
                //cout<<Row_Index<<" "<<Col_Index<<" "<<endl;
                PFT_Type_Index = int(PFT_Dominant_Index(Row_Index,Col_Index));
                //cout<<PFT_Type_Index<<endl;
                //cout<<Parameter_PFT_Space_Ensemble_SubBlock_Sens(PFT_Type_Index,Row_Index,Col_Index)<<endl;
                Parameter_Space_SubBlock(Range::all(),Row_Index,Col_Index) = Parameter_PFT_Space_Ensemble_SubBlock_Sens(Range::all(),PFT_Type_Index,Row_Index,Col_Index);
            }
        }
    """
    
    args =['-pthread -O3 -fopenmp']
    if socket.gethostname()[0:4] == 'node':
        args =['-pthread -O3 -fopenmp']
    if os.name == 'nt' or os.name == 'posix':
        libs=['gomp']
    if os.name == 'posix' and socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if os.name == 'posix' and socket.gethostname()[0:4] == 'node':
        libs=['gomp','blas','lapack']
    dirs1=['.',DAS_Depends_Path + 'include',DAS_Depends_Path + 'include/python2.7']
    dirs2=['.',DAS_Depends_Path + 'lib',DAS_Depends_Path + 'lib64']
    vars_list = ['Row_Numbers','Col_Numbers','maxpft','Parameter_Space_SubBlock','PFT_Dominant_Index','Parameter_PFT_Space_Ensemble_SubBlock_Sens','omp_get_num_procs_ParFor']
    scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])
    if socket.gethostname()[0:4] == 'node':
        scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])

    #print numpy.min(Parameter_Space_SubBlock),numpy.max(Parameter_Space_SubBlock)
    
    return

def ParFor_Sand_Clay_Organic(Sand_Top_Region,Sand_Sub_Region,Clay_Top_Region,Clay_Sub_Region,\
                         Organic_Top_Region,Organic_Sub_Region,
                         Row_Numbers,Col_Numbers,DAS_Depends_Path,omp_get_num_procs_ParFor):
    
    Sand_Top_Region = numpy.asarray(Sand_Top_Region,dtype=numpy.float32)
    Sand_Sub_Region = numpy.asarray(Sand_Sub_Region,dtype=numpy.float32)
    Clay_Top_Region = numpy.asarray(Clay_Top_Region,dtype=numpy.float32)
    Clay_Sub_Region = numpy.asarray(Clay_Sub_Region,dtype=numpy.float32)
    Organic_Top_Region = numpy.asarray(Organic_Top_Region,dtype=numpy.float32)
    Organic_Sub_Region = numpy.asarray(Organic_Sub_Region,dtype=numpy.float32)
    
    Sand_Mat = numpy.zeros((10,Row_Numbers,Col_Numbers),dtype=numpy.float32)
    Clay_Mat = numpy.zeros((10,Row_Numbers,Col_Numbers),dtype=numpy.float32)
    Organic_Mat = numpy.zeros((10,Row_Numbers,Col_Numbers),dtype=numpy.float32)
    
    support_code="""
        #define BZ_THREADSAFE
        #include <blitz/array.h>
        #include <algorithm>
    """
    code="""
        using namespace std;
        using namespace blitz;
        
        int Row_Index, Col_Index, Soil_Layer_Index;
        
        double Sand_Increment,Clay_Increment,Organic_Increment;
        
        omp_set_num_threads(omp_get_num_procs_ParFor);
        
        #pragma omp parallel for schedule(static) collapse(2) private(Row_Index,Col_Index,Soil_Layer_Index,Sand_Increment,Clay_Increment,Organic_Increment) shared(Sand_Mat,Clay_Mat,Organic_Mat)
        for(Row_Index=0;Row_Index<Row_Numbers;Row_Index++){
            for(Col_Index=0;Col_Index<Col_Numbers;Col_Index++){
                Sand_Increment = (Sand_Sub_Region(Row_Index,Col_Index) - Sand_Top_Region(Row_Index,Col_Index)) / 10.0;
                Clay_Increment = (Clay_Sub_Region(Row_Index,Col_Index) - Clay_Top_Region(Row_Index,Col_Index)) / 10.0;
                for(Soil_Layer_Index=0; Soil_Layer_Index < 10; Soil_Layer_Index++){
                    Sand_Mat(Soil_Layer_Index,Row_Index,Col_Index) = Sand_Top_Region(Row_Index,Col_Index) + Sand_Increment * Soil_Layer_Index;
                    Clay_Mat(Soil_Layer_Index,Row_Index,Col_Index) = Clay_Top_Region(Row_Index,Col_Index) + Clay_Increment * Soil_Layer_Index;
                }
                Organic_Increment = (Organic_Sub_Region(Row_Index,Col_Index) - Organic_Top_Region(Row_Index,Col_Index)) / 8.0;
                for(Soil_Layer_Index=0; Soil_Layer_Index < 8; Soil_Layer_Index++){
                    Organic_Mat(Soil_Layer_Index,Row_Index,Col_Index) = Organic_Top_Region(Row_Index,Col_Index) + Organic_Increment * Soil_Layer_Index;
                }
            }
        }
    """
    
    args =['-pthread -O3 -fopenmp']
    libs=['gomp']
    if socket.gethostname()[0:4] == 'node':
        args =['-pthread -O3 -fopenmp']
    if os.name == 'nt' or os.name == 'posix':
        libs=['gomp']
    if os.name == 'posix' and socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if os.name == 'posix' and socket.gethostname()[0:4] == 'node':
        libs=['gomp','blas','lapack']
    dirs1=['.',DAS_Depends_Path + 'include',DAS_Depends_Path + 'include/python2.7']
    dirs2=['.',DAS_Depends_Path + 'lib',DAS_Depends_Path + 'lib64']
    vars_list = ['Sand_Top_Region','Sand_Sub_Region','Clay_Top_Region','Clay_Sub_Region','Organic_Top_Region','Organic_Sub_Region', \
            'Sand_Mat','Clay_Mat','Organic_Mat','Row_Numbers','Col_Numbers','DAS_Depends_Path','omp_get_num_procs_ParFor']
    scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])
    if socket.gethostname()[0:4] == 'node':
        scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])
    
    print "Sand_Mat",numpy.min(Sand_Mat),numpy.max(Sand_Mat)
    print "Clay_Mat",numpy.min(Clay_Mat),numpy.max(Clay_Mat)
    print "Organic_Mat",numpy.min(Organic_Mat[0:8,:,:]),numpy.max(Organic_Mat[0:8,:,:])
    
    return Sand_Mat,Clay_Mat,Organic_Mat


def ParFor_H_Operator(H_Out,ny,Mask_False,Obs_Index,Soil_Moisture_DA_Flag,\
                      SensorVariable,SensorType,Soil_Layer_Index_DA,State_DIM_Single_Layer,\
                      Parameter_Optimization_Flag,Bias_Estimation_Option_Model, Bias_Estimation_Option_Obs,DAS_Depends_Path,omp_get_num_procs_ParFor):
    
    support_code="""
        #define BZ_THREADSAFE
        #include <blitz/array.h>
    """
    code="""
        
        using namespace std;
        using namespace blitz;
        
        int ny_index, H_Col_Index, Mask_Index;
                
        omp_set_num_threads(omp_get_num_procs_ParFor);
        omp_set_num_threads(1);
        
        Array<TinyVector<int,1>,1> Find_Vec;
        
        #pragma omp parallel for schedule(static) private(Mask_Index, ny_index, H_Col_Index, Find_Vec) shared(H_Out)  
        for(ny_index=0;ny_index<ny;ny_index++){
            find(Find_Vec,Mask_False == Obs_Index(ny_index));
            //cout<<"ny_index"<<" "<<ny_index<<" "<<Find_Vec<<"Find_Vec(0) "<<Find_Vec(0)<<"Find_Vec(0)(0) "<<Find_Vec(0)(0)<<endl;
            
            if(Find_Vec.size() > 0){
                if((!Parameter_Optimization_Flag) && (!(Bias_Estimation_Option_Model == 1 || Bias_Estimation_Option_Obs == 1))){
                    if((Soil_Moisture_DA_Flag == 1) && (SensorVariable == "Soil_Moisture")){
                        if(SensorType == "InSitu"){
                            H_Col_Index = Find_Vec(0)(0) + Soil_Layer_Index_DA*State_DIM_Single_Layer;
                        }
                        else{
                            H_Col_Index = Find_Vec(0)(0);
                        }
                    }
                    else{
                        H_Col_Index = Find_Vec(0)(0);
                    }
                }
                else if((!Parameter_Optimization_Flag) && (Bias_Estimation_Option_Model == 1 || Bias_Estimation_Option_Obs == 1)){
                    H_Col_Index = Find_Vec(0)(0);
                }
                else{
                    H_Col_Index = Find_Vec(0)(0);
                }
                Find_Vec.free();
                //cout<<"H_Col_Index "<<H_Col_Index<<" "<<endl;
                H_Out(ny_index,H_Col_Index) = 1;
            }
        }
    """
    
    args =['-pthread -O3 -fopenmp']
    if socket.gethostname()[0:4] == 'node':
        args =['-pthread -O3 -fopenmp']
    if os.name == 'nt' or os.name == 'posix':
        libs=['gomp']
    if os.name == 'posix' and socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if os.name == 'posix' and socket.gethostname()[0:4] == 'node':
        libs=['gomp','blas','lapack']
    dirs1=['.',DAS_Depends_Path + 'include',DAS_Depends_Path + 'include/python2.7']
    dirs2=['.',DAS_Depends_Path + 'lib',DAS_Depends_Path + 'lib64']
    vars_list = ['ny','Mask_False','Obs_Index','Soil_Moisture_DA_Flag','SensorVariable','SensorType',\
                 'Soil_Layer_Index_DA','State_DIM_Single_Layer','Parameter_Optimization_Flag','Bias_Estimation_Option_Model','Bias_Estimation_Option_Obs','omp_get_num_procs_ParFor','H_Out']
    scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])
    if socket.gethostname()[0:4] == 'node':
        scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])

def ParFor_Texture_Check(Dim_Soil_Par, Ensemble_Number, Row_Numbers, Col_Numbers, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, Parameter_Soil_Space_Ensemble, DAS_Depends_Path, omp_get_num_procs_ParFor):
    
    Parameter_Soil_Space_Ensemble_Out = numpy.asarray(Parameter_Soil_Space_Ensemble)
    
    support_code="""
        #define BZ_THREADSAFE
        #include <blitz/array.h>
    """
    code="""
        using namespace std;
        using namespace blitz;
        
        int Ens_Index, Soil_Layer_Index_Sub,  Row_Index, Col_Index;
        double Texture_Sum, Ratio, Diff, Diff_Part1, Diff_Part2;
        
        omp_set_num_threads(omp_get_num_procs_ParFor);
        
        #pragma omp parallel private(Ens_Index, Soil_Layer_Index_Sub, Row_Index, Col_Index, Texture_Sum, Ratio, Diff, Diff_Part1, Diff_Part2) shared(Parameter_Soil_Space_Ensemble_Out)
        {
            #pragma omp for schedule(static)
            for(Ens_Index=0;Ens_Index<Ensemble_Number;Ens_Index++){
                cout<<flush;
                for(Row_Index=0;Row_Index<Row_Numbers;Row_Index++){
                    for(Col_Index=0;Col_Index<Col_Numbers;Col_Index++){
                        for(Soil_Layer_Index_Sub=0;Soil_Layer_Index_Sub<Soil_Texture_Layer_Opt_Num;Soil_Layer_Index_Sub++){
                            Texture_Sum = (Parameter_Soil_Space_Ensemble_Out(Ens_Index,Soil_Layer_Index_Sub,Row_Index,Col_Index) + Parameter_Soil_Space_Ensemble_Out(Ens_Index,Soil_Layer_Index_Sub+Soil_Texture_Layer_Opt_Num,Row_Index,Col_Index));
                            if(Texture_Sum > 98.0){
                                Ratio = Parameter_Soil_Space_Ensemble_Out(Ens_Index,Soil_Layer_Index_Sub,Row_Index,Col_Index) / (Parameter_Soil_Space_Ensemble_Out(Ens_Index,Soil_Layer_Index_Sub,Row_Index,Col_Index)+Parameter_Soil_Space_Ensemble_Out(Ens_Index,Soil_Layer_Index_Sub+Soil_Texture_Layer_Opt_Num,Row_Index,Col_Index));
                                Diff = Texture_Sum - 98.0;
                                Diff_Part1 = Ratio*Diff;
                                Diff_Part2 = Diff - Diff_Part1;
                                Parameter_Soil_Space_Ensemble_Out(Ens_Index,Soil_Layer_Index_Sub,Row_Index,Col_Index) -= Diff_Part1;
                                Parameter_Soil_Space_Ensemble_Out(Ens_Index,Soil_Layer_Index_Sub+Soil_Texture_Layer_Opt_Num,Row_Index,Col_Index) -= Diff_Part2;
                            }
                            //if(Texture_Sum < Soil_Sand_Clay_Sum(Soil_Layer_Index_Sub,Row_Index,Col_Index)){
                            //    Ratio = Parameter_Soil_Space_Ensemble_Out(Ens_Index,Soil_Layer_Index_Sub,Row_Index,Col_Index) / (Parameter_Soil_Space_Ensemble_Out(Ens_Index,Soil_Layer_Index_Sub,Row_Index,Col_Index)+Parameter_Soil_Space_Ensemble(Ens_Index,Soil_Layer_Index_Sub+Soil_Texture_Layer_Opt_Num,Row_Index,Col_Index));
                            //    Diff = Soil_Sand_Clay_Sum(Soil_Layer_Index_Sub,Row_Index,Col_Index) - Texture_Sum;
                            //    Diff_Part1 = Ratio*Diff;
                            //    Diff_Part2 = Diff - Diff_Part1;
                            //    Parameter_Soil_Space_Ensemble_Out(Ens_Index,Soil_Layer_Index_Sub,Row_Index,Col_Index) += Diff_Part1;
                            //    Parameter_Soil_Space_Ensemble_Out(Ens_Index,Soil_Layer_Index_Sub+Soil_Texture_Layer_Opt_Num,Row_Index,Col_Index) += Diff_Part2;
                            //}
                        }
                    }
                }
            }
        }
    """
    
    args =['-pthread -O3 -fopenmp']
    libs=['gomp']
    if socket.gethostname()[0:4] == 'node':
        args =['-pthread -O3 -fopenmp']
    if os.name == 'nt' or os.name == 'posix':
        libs=['gomp']
    if os.name == 'posix' and socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if os.name == 'posix' and socket.gethostname()[0:4] == 'node':
        libs=['gomp','blas','lapack']
    dirs1=['.',DAS_Depends_Path + 'include',DAS_Depends_Path + 'include/python2.7']
    dirs2=['.',DAS_Depends_Path + 'lib',DAS_Depends_Path + 'lib64']
    vars_list = ['Ensemble_Number','Soil_Texture_Layer_Opt_Num','Soil_Sand_Clay_Sum','Row_Numbers','Col_Numbers','Parameter_Soil_Space_Ensemble','Parameter_Soil_Space_Ensemble_Out','omp_get_num_procs_ParFor']
    scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])
    if socket.gethostname()[0:4] == 'node':
        scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])

    return Parameter_Soil_Space_Ensemble_Out

def ParFor_Ratio(CLM_NA,ColumnNum,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,Array_Origin,ratio_upper,ratio_lower,Mask_Region,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue):
    
    varout = numpy.zeros(ColumnNum,dtype=numpy.float32)
    
    support_code="""
        #define BZ_THREADSAFE
        #include <blitz/array.h>
    """
    code="""
        using namespace std;
        using namespace blitz;
        
        int Column_Index,  Row_Index, Col_Index;
        
        omp_set_num_threads(omp_get_num_procs_ParFor);
        
        #pragma omp parallel private(Column_Index,Row_Index,Col_Index) shared(varout)
        {
            #pragma omp for schedule(static)
            for(Column_Index=0;Column_Index<ColumnNum;Column_Index++){
                cout<<flush;
                Row_Index=Row_Numbers-cols1d_jxy(Column_Index);
                Col_Index=cols1d_ixy(Column_Index) - 1;
                varout(Column_Index)=Array_Origin(Column_Index);
                if(!Mask_Region(Row_Index,Col_Index)){
                    varout(Column_Index) = Array_Origin(Column_Index) * ratio_upper(Column_Index) / ratio_lower(Column_Index);
                }
            }
        }
    """
    
    args =['-pthread -O3 -fopenmp']
    if socket.gethostname()[0:4] == 'node':
        args =['-pthread -O3 -fopenmp']
    if os.name == 'nt' or os.name == 'posix':
        libs=['gomp']
    if os.name == 'posix' and socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if os.name == 'posix' and socket.gethostname()[0:4] == 'node':
        libs=['gomp','blas','lapack']
    dirs1=['.',DAS_Depends_Path + 'include',DAS_Depends_Path + 'include/python2.7']
    dirs2=['.',DAS_Depends_Path + 'lib',DAS_Depends_Path + 'lib64']
    vars_list = ['CLM_NA','ColumnNum','Row_Numbers','Col_Numbers','cols1d_ixy','cols1d_jxy','Array_Origin','ratio_upper','ratio_lower','Mask_Region','varout','omp_get_num_procs_ParFor','NAvalue']
    scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])
    if socket.gethostname()[0:4] == 'node':
        scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])

    #print numpy.min(varout),numpy.max(varout)
    return varout


def ParFor_Soil_Moisture(CLM_NA,ColumnNum,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,Soil_Layer_Index_DA,Snow_Layer_Num,Soil_Layer_Num,CLM_Soil_Moisture,CLM_Soil_Layer_Thickness,Density_of_liquid_water,Density_of_ice,\
                         Teta_Saturated,cols1d_ityplun,CLM_Soil_Temperature,Freezing_temperature_of_fresh_water,H2OSOI_LIQ_In,Mask_Region,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue):
    
    H2OSOI_LIQ = numpy.zeros_like(H2OSOI_LIQ_In)
    H2OSOI_LIQ[::] = H2OSOI_LIQ_In[::]
    
    support_code="""
        #define BZ_THREADSAFE
        #include <blitz/array.h>
        #include <algorithm>
    """
    code="""
        using namespace std;
        using namespace blitz;
        
        int Column_Index,  Row_Index, Col_Index, Soil_Layer_Index;
        int cols1d_ityplun_index;
        double Diff_Temp1, Diff_Temp2;
        
        omp_set_num_threads(omp_get_num_procs_ParFor);
        
        #pragma omp parallel private(Column_Index,Row_Index,Col_Index, cols1d_ityplun_index, Soil_Layer_Index, Diff_Temp1, Diff_Temp2) shared(H2OSOI_LIQ, CLM_Soil_Moisture,CLM_Soil_Layer_Thickness,Teta_Saturated)
        {
            #pragma omp for schedule(static)
            for(Column_Index=0;Column_Index<ColumnNum;Column_Index++){
                cout<<flush;
                Row_Index=Row_Numbers-cols1d_jxy(Column_Index);
                Col_Index=cols1d_ixy(Column_Index) - 1;
                cols1d_ityplun_index = cols1d_ityplun(Column_Index);
                
                if(cols1d_ityplun_index != 3){
                
                    if(!Mask_Region(Row_Index,Col_Index)){
                        for(Soil_Layer_Index=0; Soil_Layer_Index < Soil_Layer_Num; Soil_Layer_Index++){
                            H2OSOI_LIQ(Column_Index,Snow_Layer_Num+Soil_Layer_Index) = H2OSOI_LIQ_In(Column_Index,Snow_Layer_Num+Soil_Layer_Index);
                            //H2OSOI_ICE(Column_Index,Snow_Layer_Num+Soil_Layer_Index) = H2OSOI_ICE_In(Column_Index,Snow_Layer_Num+Soil_Layer_Index);
                        }
                        
                        if(cols1d_ityplun_index == 1 || cols1d_ityplun_index == 5 || cols1d_ityplun_index == 8){
                            // Transform the CLM Soil Moisture Unit (m3/m3) to (kg/m2) after data assimilation
                            for(Soil_Layer_Index=0; Soil_Layer_Index < Soil_Layer_Num; Soil_Layer_Index++){
                                if(CLM_Soil_Moisture(Soil_Layer_Index,Row_Index,Col_Index) != NAvalue && CLM_Soil_Layer_Thickness(Soil_Layer_Index,Row_Index,Col_Index) < CLM_NA){
                                    if(CLM_Soil_Temperature(Soil_Layer_Index,Row_Index, Col_Index) <= Freezing_temperature_of_fresh_water){
                                        //H2OSOI_ICE(Column_Index,Snow_Layer_Num+Soil_Layer_Index) = CLM_Soil_Moisture(Soil_Layer_Index,Row_Index,Col_Index) * (CLM_Soil_Layer_Thickness(Soil_Layer_Index,Row_Index,Col_Index) * double(Density_of_ice));
                                        //H2OSOI_LIQ(Column_Index,Snow_Layer_Num+Soil_Layer_Index) = 0.0;
                                        H2OSOI_LIQ(Column_Index,Snow_Layer_Num+Soil_Layer_Index) = CLM_Soil_Moisture(Soil_Layer_Index,Row_Index,Col_Index) * (CLM_Soil_Layer_Thickness(Soil_Layer_Index,Row_Index,Col_Index) * double(Density_of_liquid_water));
                                        // Update the Soil Ice
                                        //if(CLM_Soil_Moisture(Soil_Layer_Index,Row_Index,Col_Index) + CLM_Soil_Ice(Soil_Layer_Index,Row_Index,Col_Index) > Teta_Saturated(Row_Index,Col_Index)){
                                           //Diff_Temp1 = Teta_Saturated(Row_Index,Col_Index) - CLM_Soil_Moisture(Soil_Layer_Index,Row_Index,Col_Index);
                                            //CLM_Soil_Ice(Snow_Layer_Num+Soil_Layer_Index,Row_Index,Col_Index) = max(Diff_Temp1,0.0);
                                            //H2OSOI_ICE(Column_Index,Snow_Layer_Num+Soil_Layer_Index) = CLM_Soil_Ice(Soil_Layer_Index,Row_Index,Col_Index) * (CLM_Soil_Layer_Thickness(Soil_Layer_Index,Row_Index,Col_Index) * double(Density_of_ice));
                                        //}
                                        //for(Soil_Layer_Index=0; Soil_Layer_Index < Soil_Layer_Num; Soil_Layer_Index++){
                                        //    H2OSOI_LIQ(Column_Index,Snow_Layer_Num+Soil_Layer_Index) = H2OSOI_LIQ_In(Column_Index,Snow_Layer_Num+Soil_Layer_Index) * H2OSOI_LIQ(Column_Index,Snow_Layer_Num+Soil_Layer_Index_DA) / H2OSOI_LIQ_In(Column_Index,Snow_Layer_Num+Soil_Layer_Index_DA);
                                        //    H2OSOI_ICE(Column_Index,Snow_Layer_Num+Soil_Layer_Index) = H2OSOI_ICE_In(Column_Index,Snow_Layer_Num+Soil_Layer_Index) * H2OSOI_ICE(Column_Index,Snow_Layer_Num+Soil_Layer_Index_DA) / H2OSOI_ICE_In(Column_Index,Snow_Layer_Num+Soil_Layer_Index_DA);
                                        //}
                                    }
                                    else{
                                        H2OSOI_LIQ(Column_Index,Snow_Layer_Num+Soil_Layer_Index) = CLM_Soil_Moisture(Soil_Layer_Index,Row_Index,Col_Index) * (CLM_Soil_Layer_Thickness(Soil_Layer_Index,Row_Index,Col_Index) * double(Density_of_liquid_water));
                                        //H2OSOI_ICE(Column_Index,Snow_Layer_Num+Soil_Layer_Index) = 0.0;
                                        //for(Soil_Layer_Index=0; Soil_Layer_Index < Soil_Layer_Num; Soil_Layer_Index++){
                                        //    H2OSOI_LIQ(Column_Index,Snow_Layer_Num+Soil_Layer_Index) = H2OSOI_LIQ_In(Column_Index,Snow_Layer_Num+Soil_Layer_Index) * H2OSOI_LIQ(Column_Index,Snow_Layer_Num+Soil_Layer_Index_DA) / H2OSOI_LIQ_In(Column_Index,Snow_Layer_Num+Soil_Layer_Index_DA);
                                        //    if(CLM_Soil_Temperature(Soil_Layer_Index,Row_Index, Col_Index) > Freezing_temperature_of_fresh_water){
                                        //        H2OSOI_ICE(Column_Index,Snow_Layer_Num+Soil_Layer_Index) = 0.0;
                                        //    }
                                        //}
                                    }
                                }
                            }
                        }
                        else
                        {
                            continue;
                        }
                    }
                }
            }
        }
    """
    
    args =['-pthread -O3 -fopenmp']
    if socket.gethostname()[0:4] == 'node':
        args =['-pthread -O3 -fopenmp']
    if os.name == 'nt' or os.name == 'posix':
        libs=['gomp']
    if os.name == 'posix' and socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if os.name == 'posix' and socket.gethostname()[0:4] == 'node':
        libs=['gomp','blas','lapack']
    dirs1=['.',DAS_Depends_Path + 'include',DAS_Depends_Path + 'include/python2.7']
    dirs2=['.',DAS_Depends_Path + 'lib',DAS_Depends_Path + 'lib64']
    vars_list = ['CLM_NA','ColumnNum','Row_Numbers','Col_Numbers','cols1d_ixy','cols1d_jxy','Soil_Layer_Index_DA','Snow_Layer_Num','Soil_Layer_Num','CLM_Soil_Moisture','CLM_Soil_Layer_Thickness',\
            'Density_of_liquid_water','Density_of_ice','Teta_Saturated','cols1d_ityplun','CLM_Soil_Temperature',\
            'Freezing_temperature_of_fresh_water','H2OSOI_LIQ_In','H2OSOI_LIQ','Mask_Region','omp_get_num_procs_ParFor','NAvalue']
    scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])
    if socket.gethostname()[0:4] == 'node':
        scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])

    #print numpy.shape(H2OSOI_LIQ[:,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)]), numpy.shape(H2OSOI_ICE[:,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)])
    
    return H2OSOI_LIQ[:,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)]

def ParFor_Check_Outliers(Land_Mask_Data, Matrix_In,Row_Numbers,Col_Numbers,Matrix_Min,Matrix_Max,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue):
    
    support_code="""
        #define BZ_THREADSAFE
        #include <blitz/array.h>
    """
    code="""
        using namespace std;
        using namespace blitz;
        
        int Row_Index, Col_Index;
        
        omp_set_num_threads(omp_get_num_procs_ParFor);
        
        #pragma omp parallel private(Row_Index,Col_Index) shared(Matrix_In)
        {
            #pragma omp for schedule(static) collapse(2)
            for(Row_Index=0;Row_Index<Row_Numbers;Row_Index++){
                for(Col_Index=0;Col_Index<Col_Numbers;Col_Index++){
                    cout<<flush;
                    if(Land_Mask_Data(Row_Index,Col_Index) != NAvalue){
                        Matrix_In(Row_Index,Col_Index) = Matrix_In(Row_Index,Col_Index);
                        if(Matrix_In(Row_Index,Col_Index) < Matrix_Min(Row_Index,Col_Index)){
                            Matrix_In(Row_Index,Col_Index) = Matrix_Min(Row_Index,Col_Index);
                        }
                        else if(Matrix_In(Row_Index,Col_Index) > Matrix_Max(Row_Index,Col_Index)){
                            //cout<<Matrix_In(Row_Index,Col_Index)<<endl;
                            Matrix_In(Row_Index,Col_Index) = Matrix_Max(Row_Index,Col_Index);
                        }
                    }
                    else{
                        Matrix_In(Row_Index,Col_Index) = NAvalue;
                    }
                }
            }
        }
        
    """
    
    args =['-pthread -O3 -fopenmp']
    if socket.gethostname()[0:4] == 'node':
        args =['-pthread -O3 -fopenmp']
    if os.name == 'nt' or os.name == 'posix':
        libs=['gomp']
    if os.name == 'posix' and socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if os.name == 'posix' and socket.gethostname()[0:4] == 'node':
        libs=['gomp','blas','lapack']
    dirs1=['.',DAS_Depends_Path + 'include',DAS_Depends_Path + 'include/python2.7']
    dirs2=['.',DAS_Depends_Path + 'lib',DAS_Depends_Path + 'lib64']
    vars_list = ['Land_Mask_Data','Matrix_In','Row_Numbers','Col_Numbers','Matrix_Min','Matrix_Max','DAS_Depends_Path','omp_get_num_procs_ParFor','NAvalue']
    scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])
    if socket.gethostname()[0:4] == 'node':
        scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])

    return

def ParFor_Check_Outliers_NA(Land_Mask_Data, Matrix_In,Row_Numbers,Col_Numbers,Matrix_Min,Matrix_Max,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue):
    
    support_code="""
        #define BZ_THREADSAFE
        #include <blitz/array.h>
    """
    code="""
        using namespace std;
        using namespace blitz;
        
        int Row_Index, Col_Index;
        
        omp_set_num_threads(omp_get_num_procs_ParFor);
        
        #pragma omp parallel private(Row_Index,Col_Index) shared(Matrix_In)
        {
            #pragma omp for schedule(static) collapse(2)
            
            for(Row_Index=0;Row_Index<Row_Numbers;Row_Index++){
                for(Col_Index=0;Col_Index<Col_Numbers;Col_Index++){
                    cout<<flush;
                    if(Land_Mask_Data(Row_Index,Col_Index) != NAvalue){
                        Matrix_In(Row_Index,Col_Index) = Matrix_In(Row_Index,Col_Index);
                        if(Matrix_In(Row_Index,Col_Index) < Matrix_Min(Row_Index,Col_Index)){
                            Matrix_In(Row_Index,Col_Index) = NAvalue;
                        }
                        else if(Matrix_In(Row_Index,Col_Index) > Matrix_Max(Row_Index,Col_Index)){
                            Matrix_In(Row_Index,Col_Index) = NAvalue;
                        }
                    }
                    else{
                        Matrix_In(Row_Index,Col_Index) = NAvalue;
                    }
                    //cout<<"Row_Index"<<Row_Index<<"Col_Index"<<Col_Index<<endl;
                    //cout<<Matrix_In(Row_Index,Col_Index)<<endl;
                }
            }
        }
        
    """
    
    args =['-pthread -O3 -fopenmp']
    if socket.gethostname()[0:4] == 'node':
        args =['-pthread -O3 -fopenmp']
    if os.name == 'nt' or os.name == 'posix':
        libs=['gomp']
    if os.name == 'posix' and socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if socket.gethostname()[0] == 'j':
        libs=['gomp','blas','lapack']
    if os.name == 'posix' and socket.gethostname()[0:4] == 'node':
        libs=['gomp','blas','lapack']
    dirs1=['.',DAS_Depends_Path + 'include',DAS_Depends_Path + 'include/python2.7']
    dirs2=['.',DAS_Depends_Path + 'lib',DAS_Depends_Path + 'lib64']
    vars_list = ['Land_Mask_Data','Matrix_In','Row_Numbers','Col_Numbers','Matrix_Min','Matrix_Max','DAS_Depends_Path','omp_get_num_procs_ParFor','NAvalue']
    scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])
    if socket.gethostname()[0:4] == 'node':
        scipy.weave.inline(code,arg_names = vars_list,support_code=support_code,extra_compile_args=args,include_dirs=dirs1,libraries=libs,library_dirs=dirs2, type_converters=scipy.weave.converters.blitz,compiler='gcc',verbose=2,headers=['<omp.h>'])

    return

