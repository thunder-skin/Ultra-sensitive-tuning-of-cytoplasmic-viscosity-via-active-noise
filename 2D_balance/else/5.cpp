#include <iostream>
#include <random>
#include <chrono>
#include <ctime>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

//定义x和y方向格数，总步数计数
int x_num=0;
int y_num=0;
int count=0;
std::string filename="1500";
int particle_num=16744;

//定义物理界参数
double x_sep=0;
double y_sep=0;
int multi_time=3;
int x_tot=100*multi_time;
int y_tot=100*multi_time;
float k=1;
float eta=1;
double sheer=0;
double t_prim=1e-2;
double t_max=1e-1;
double s_prim=1e-3*x_tot;
double amp_control=1e-3;

//定义工作
bool prdiagnose=false;
std::vector<double> force_new(2*particle_num,0);
std::vector<double> force_old(2*particle_num,0);
std::vector<double> velocy_new(2*particle_num,0);
std::vector<double> velocy_old(2*particle_num,0);
std::vector<double> place_data(2*particle_num,0);




//定义物理世界
std::vector<float> radius={1,1.4};
float max_r=1.4;
int balance_max=1e4;
int random_seed=std::stoi(filename)+10000;
std::mt19937 gen(random_seed);
std::random_device rd;

//粒子类
class Particle{
public:
    float r;
    double x;
    double y;
    int id;
    int x_cell;
    int y_cell;
    Particle():r(0),x(0),y(0),id(0){}
    Particle(float radius,double particle_x,double particle_y,int particle_id)
        :r(radius),x(particle_x),y(particle_y),id(particle_id){
    x_cell=x/x_sep;
    y_cell=y/y_sep;
    }
};

//创造随机初始世界
std::vector<std::vector<std::vector<Particle>>> generate_world() {
    double cell=2*max_r+0.1;
    x_num=x_tot/cell;
    y_num=y_tot/cell;
    x_sep=x_tot/double(x_num);
    y_sep=y_tot/double(y_num);
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    std::vector<std::vector<std::vector<Particle>>> savings(x_num, std::vector<std::vector<Particle>>(y_num, std::vector<Particle>()));
    for (int i=0;i<particle_num;++i) {
        int choose_id=i%2;
        float r_i=radius[choose_id];
        double random_num=0;
        random_num=dis(gen);
        double x_i=x_tot*random_num;
        random_num=dis(gen);
        double y_i=y_tot*random_num;
        Particle particle_i(r_i,x_i,y_i,i);
        savings[particle_i.x_cell][particle_i.y_cell].push_back(particle_i);
    }
    return savings;
}

//实例化随机初始世界
std::vector<std::vector<std::vector<Particle>>> savings=generate_world();

//计算一对粒子受力
std::vector<double> calculate_force(double x_self,double y_self,double x,double y,float r_self,float r,int id1,int id2){
    if (x_self-x>2*x_tot/3){
        x_self-=x_tot;
    }
    else if (x-x_self>2*x_tot/3){
        x_self+=x_tot;
    }
    if ((y_self>2*y_tot/3) and (y<y_tot/3)){
        y_self-=y_tot;
    }
    else if ((y_self<y_tot/3) and (y>2*y_tot/3)){
        y_self+=y_tot;
    }
    if ((pow((x_self-x),2)+pow((y_self-y),2)>pow((r_self+r),2)) or id1==id2) {
        std::vector<double> force(3, 0);
        return force;
    }
    else{
        double delta_r=pow(pow((x_self-x),2)+pow((y_self-y),2),0.5);
        double f_tot=k*(r_self+r-delta_r);
        double fx=f_tot*(x_self-x)/delta_r;
        double fy=f_tot*(y_self-y)/delta_r;
        double Lx=fx*(y-y_self);
        std::vector<double> force{fx,fy,Lx};
        return force;
    }


}

//计算一对粒子能量
double calculate_energy(double x_self,double y_self,double x,double y,float r_self,float r,int id1,int id2){
    if (x_self-x>2*x_tot/3){
        x_self-=x_tot;
    }
    else if (x-x_self>2*x_tot/3){
        x_self+=x_tot;
    }
    if ((y_self>2*y_tot/3) and (y<y_tot/3)){
        y_self-=y_tot;
    }
    else if ((y_self<y_tot/3) and (y>2*y_tot/3)){
        y_self+=y_tot;
    }
    if ((pow((x_self-x),2)+pow((y_self-y),2)>pow((r_self+r),2)) or id1==id2) {
        return 0;
    }
    else{
        double delta_r=pow(pow((x_self-x),2)+pow((y_self-y),2),0.5);
        double E_tot=k*(r_self+r-delta_r)*(r_self+r-delta_r)*0.5;
        return E_tot;
    }
}

//寻找临近受力
std::vector<double> neighbor_force(Particle part_self){
    double force_x=0;
    double force_y=0;
    double L_tot=0;
    int f_cell=sheer/x_sep;
    for (int i=-1;i<=1;++i){
        for (int j=0;j<savings[fmod((part_self.x_cell+i+x_num),x_num)][part_self.y_cell].size();++j){
            Particle part_else=savings[fmod((part_self.x_cell+i+x_num),x_num)][part_self.y_cell][j];
            std::vector<double> single_force=calculate_force(part_self.x,part_self.y,part_else.x,part_else.y,part_self.r,part_else.r,part_self.id,part_else.id);
            force_x+=single_force[0];
            force_y+=single_force[1];
            L_tot+=single_force[2];
        }
    }
    if (part_self.y_cell==0){
        for (int i=-1;i<=1;++i){
            for (int j=0;j<savings[fmod((part_self.x_cell+i+x_num),x_num)][1].size();++j){
                Particle part_else=savings[fmod((part_self.x_cell+i+x_num),x_num)][1][j];
                std::vector<double> single_force=calculate_force(part_self.x,part_self.y,part_else.x,part_else.y,part_self.r,part_else.r,part_self.id,part_else.id);
                force_x+=single_force[0];
                force_y+=single_force[1];
                L_tot+=single_force[2];
            }
        }
        for (int i=-2;i<=2;++i){
            for (int j=0;j<savings[fmod((part_self.x_cell+i+x_num+f_cell),x_num)][y_num-1].size();++j){
                Particle part_else=savings[fmod((part_self.x_cell+i+x_num+f_cell),x_num)][y_num-1][j];
                std::vector<double> single_force=calculate_force(part_self.x,part_self.y,part_else.x-sheer,part_else.y,part_self.r,part_else.r,part_self.id,part_else.id);
                force_x+=single_force[0];
                force_y+=single_force[1];
                L_tot+=single_force[2];
            }
        }
    }
    else if (part_self.y_cell==y_num-1){
        for (int i=-1;i<=1;++i){
            for (int j=0;j<savings[fmod((part_self.x_cell+i+x_num),x_num)][y_num-2].size();++j){
                Particle part_else=savings[fmod((part_self.x_cell+i+x_num),x_num)][y_num-2][j];
                std::vector<double> single_force=calculate_force(part_self.x,part_self.y,part_else.x,part_else.y,part_self.r,part_else.r,part_self.id,part_else.id);
                force_x+=single_force[0];
                force_y+=single_force[1];
                L_tot+=single_force[2];
            }
        }
        for (int i=-2;i<=2;++i){
            for (int j=0;j<savings[fmod((part_self.x_cell+i+x_num-f_cell),x_num)][0].size();++j){
                Particle part_else=savings[fmod((part_self.x_cell+i+x_num-f_cell),x_num)][0][j];
                std::vector<double> single_force=calculate_force(part_self.x,part_self.y,part_else.x+sheer,part_else.y,part_self.r,part_else.r,part_self.id,part_else.id);
                force_x+=single_force[0];
                force_y+=single_force[1];
                L_tot+=single_force[2];
            }
        }
    }
    else{
        for (int i=-1;i<=1;++i){
            for (int j=0;j<savings[fmod((part_self.x_cell+i+x_num),x_num)][part_self.y_cell+1].size();++j){
                Particle part_else=savings[fmod((part_self.x_cell+i+x_num),x_num)][part_self.y_cell+1][j];
                std::vector<double> single_force=calculate_force(part_self.x,part_self.y,part_else.x,part_else.y,part_self.r,part_else.r,part_self.id,part_else.id);
                force_x+=single_force[0];
                force_y+=single_force[1];
                L_tot+=single_force[2];
            }
        }
        for (int i=-1;i<=1;++i){
            for (int j=0;j<savings[fmod((part_self.x_cell+i+x_num),x_num)][part_self.y_cell-1].size();++j){
                Particle part_else=savings[fmod((part_self.x_cell+i+x_num),x_num)][part_self.y_cell-1][j];
                std::vector<double> single_force=calculate_force(part_self.x,part_self.y,part_else.x,part_else.y,part_self.r,part_else.r,part_self.id,part_else.id);
                force_x+=single_force[0];
                force_y+=single_force[1];
                L_tot+=single_force[2]; 
            }
        }
    }
    if (prdiagnose==true){
        force_new[2*part_self.id]=force_x;
        force_new[2*part_self.id+1]=force_y;
    }
    std::vector<double> force{force_x,force_y,L_tot};
    return force;
}

//计算临近能量
double neighbor_energy(Particle part_self){
    double E_all=0;
    int f_cell=sheer/x_sep;
    for (int i=-1;i<=1;++i){
        for (int j=0;j<savings[fmod((part_self.x_cell+i+x_num),x_num)][part_self.y_cell].size();++j){
            Particle part_else=savings[fmod((part_self.x_cell+i+x_num),x_num)][part_self.y_cell][j];
            double single_E=calculate_energy(part_self.x,part_self.y,part_else.x,part_else.y,part_self.r,part_else.r,part_self.id,part_else.id);
            E_all+=single_E;
        }
    }
    if (part_self.y_cell==0){
        for (int i=-1;i<=1;++i){
            for (int j=0;j<savings[fmod((part_self.x_cell+i+x_num),x_num)][1].size();++j){
                Particle part_else=savings[fmod((part_self.x_cell+i+x_num),x_num)][1][j];
                double single_E=calculate_energy(part_self.x,part_self.y,part_else.x,part_else.y,part_self.r,part_else.r,part_self.id,part_else.id);
                E_all+=single_E;
            }
        }
        for (int i=-2;i<=2;++i){
            for (int j=0;j<savings[fmod((part_self.x_cell+i+x_num+f_cell),x_num)][y_num-1].size();++j){
                Particle part_else=savings[fmod((part_self.x_cell+i+x_num+f_cell),x_num)][y_num-1][j];
                double single_E=calculate_energy(part_self.x,part_self.y,part_else.x-sheer,part_else.y,part_self.r,part_else.r,part_self.id,part_else.id);
                E_all+=single_E;
            }
        }
    }
    else if (part_self.y_cell==y_num-1){
        for (int i=-1;i<=1;++i){
            for (int j=0;j<savings[fmod((part_self.x_cell+i+x_num),x_num)][y_num-2].size();++j){
                Particle part_else=savings[fmod((part_self.x_cell+i+x_num),x_num)][y_num-2][j];
                double single_E=calculate_energy(part_self.x,part_self.y,part_else.x,part_else.y,part_self.r,part_else.r,part_self.id,part_else.id);
                E_all+=single_E;
            }
        }
        for (int i=-2;i<=2;++i){
            for (int j=0;j<savings[fmod((part_self.x_cell+i+x_num-f_cell),x_num)][0].size();++j){
                Particle part_else=savings[fmod((part_self.x_cell+i+x_num-f_cell),x_num)][0][j];
                double single_E=calculate_energy(part_self.x,part_self.y,part_else.x+sheer,part_else.y,part_self.r,part_else.r,part_self.id,part_else.id);
                E_all+=single_E;
            }
        }
    }
    else{
        for (int i=-1;i<=1;++i){
            for (int j=0;j<savings[fmod((part_self.x_cell+i+x_num),x_num)][part_self.y_cell+1].size();++j){
                Particle part_else=savings[fmod((part_self.x_cell+i+x_num),x_num)][part_self.y_cell+1][j];
                double single_E=calculate_energy(part_self.x,part_self.y,part_else.x,part_else.y,part_self.r,part_else.r,part_self.id,part_else.id);
                E_all+=single_E;
            }
        }
        for (int i=-1;i<=1;++i){
            for (int j=0;j<savings[fmod((part_self.x_cell+i+x_num),x_num)][part_self.y_cell-1].size();++j){
                Particle part_else=savings[fmod((part_self.x_cell+i+x_num),x_num)][part_self.y_cell-1][j];
                double single_E=calculate_energy(part_self.x,part_self.y,part_else.x,part_else.y,part_self.r,part_else.r,part_self.id,part_else.id);
                E_all+=single_E;
            }
        }
    }
    return E_all;
}

//计算总能量
double total_energy(){
    double Energy_total=0;
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<savings[i][j].size();++k){
                Energy_total+=neighbor_energy(savings[i][j][k]);
            }
        }
    }
    return Energy_total;
}

//存档
void write_world(double L_0){
    std::string name=std::to_string(random_seed);
    std::ofstream outFile(filename+".txt");
    outFile<<x_tot<<" "<<y_tot<<" "<<sheer<<" "<<L_0<<std::endl;
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<savings[i][j].size();++k){
                Particle part=savings[i][j][k];
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(15); // 设置浮点数输出精度
                oss << part.x << " " << part.y << " " << part.r << " " << part.id << std::endl;
                outFile << oss.str();
            }
        }
    }
    outFile.close();
}

//计算MD移动
std::vector<double> move_to(double fx,double fy,Particle part_self,float t_step){
    double x_move=fx/eta*t_step;
    double y_move=fy/eta*t_step;
    double y_final=0;
    double x_final=0;
    if (y_move+part_self.y>y_tot){
        y_final=y_move+part_self.y-y_tot;
        x_final=x_move+part_self.x-sheer;
    }
    else if (y_move+part_self.y<0){
        y_final=y_move+part_self.y+y_tot;
        x_final=x_move+part_self.x+sheer;
    }
    else{
        y_final=y_move+part_self.y;
        x_final=x_move+part_self.x;
    }
    x_final=fmod(x_final+x_tot,x_tot);
    std::vector<double> force{x_final,y_final,x_move,y_move};
    return force;
}

//MD步进和力矩
double step_forward(float t_step){
    std::vector<std::vector<std::vector<Particle>>> new_savings(x_num, std::vector<std::vector<Particle>>(y_num, std::vector<Particle>(0)));
    double sheer_force=0;
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<savings[i][j].size();++k){
                
                Particle part_old=savings[i][j][k];
                std::vector<double> force=neighbor_force(part_old);
                sheer_force+=force[2];
                std::vector<double> move=move_to(force[0],force[1],part_old,t_step);
                Particle part_new(part_old.r,move[0],move[1],part_old.id);
                new_savings[part_new.x_cell][part_new.y_cell].push_back(part_new);
            }
        }
    }
    savings=new_savings;
    return sheer_force;
}

//剪切
void build_sheer(double sheer_add){
    std::vector<std::vector<std::vector<Particle>>> new_savings(x_num, std::vector<std::vector<Particle>>(y_num, std::vector<Particle>(0)));
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<savings[i][j].size();++k){
                Particle part_old=savings[i][j][k];
                double x_new=fmod((part_old.x+part_old.y*sheer_add/y_tot+x_tot),x_tot);
                Particle part_new(part_old.r,x_new,part_old.y,part_old.id);
                new_savings[i][j].push_back(part_new);
            }
        }
    }
    savings=new_savings;
}

//建立MD平衡
double find_md_balance(int md_max){
    int count=0;
    double L_new=0;
    while (count<=md_max){
        count+=1;
        float t_step=1e2/(count+1e3);
        L_new=step_forward(t_step);
    }
    return L_new;
}

//获取受力
void active_force_vec(){
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<savings[i][j].size();++k){
                Particle partii=savings[i][j][k];
                std::vector<double> useless=neighbor_force(partii);
            }
        }
    }
}

//计算运动方向
void calculate_direction(){
    double beta_1=0;
    double beta_2=0;
    for (int k=0;k<2*particle_num;++k){
        beta_1+=force_new[k]*(force_new[k]-force_old[k]);
        beta_2+=force_old[k]*force_old[k];
    }
    for (int k=0;k<2*particle_num;++k){
        velocy_new[k]=force_new[k]+beta_1/beta_2*velocy_old[k];
    }
    velocy_old=velocy_new;
    force_old=force_new;
}

//方向迭代
double direction_iteration(){
    double E_old=1e10;
    double E_new=1e9;
    int error_times=1;
    int direc=1;
    double t_step=t_prim;
    while ((E_new/E_old)>(1+amp_control/10) or (E_new/E_old)<(1-amp_control/10)){
        
        if (E_new>E_old){
            error_times+=error_times;
            direc*=(-1);
            t_step=t_prim;
        }
        else if (t_step<t_max){
            t_step*=2;
        }
        E_old=E_new;
        for (int k=0;k<particle_num;++k){
            place_data[2*k+1]+=t_step*velocy_new[2*k+1]*direc/error_times;
            place_data[2*k]+=t_step*velocy_new[2*k]*direc/error_times;
            if (place_data[2*k+1]>=y_tot){
                place_data[2*k+1]-=y_tot;
                place_data[2*k]-=sheer;
            }
            else if (place_data[2*k+1]<0){
                place_data[2*k+1]+=y_tot;
                place_data[2*k]+=sheer;
            }
            place_data[2*k]=fmod(place_data[2*k]+x_tot,x_tot);
        }
        std::vector<std::vector<std::vector<Particle>>> new_savings(x_num, std::vector<std::vector<Particle>>(y_num, std::vector<Particle>(0)));
        for (int k=0;k<particle_num;++k){
            int choose_id=k%2;
            float r_i=radius[choose_id];
            Particle part_new(r_i,place_data[k*2],place_data[k*2+1],k);
            new_savings[part_new.x_cell][part_new.y_cell].push_back(part_new);
        }
        savings=new_savings;
        
        E_new=total_energy();
    }
    return E_new;
}

//计算总剪切矩
double total_L(){
    double L=0;
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<savings[i][j].size();++k){
                Particle part_old=savings[i][j][k];
                std::vector<double> force=neighbor_force(part_old);
                L+=force[2];
            }
        }
    }
    return L;
}

//剪切迭代
std::vector<double> sheer_itration(){
    double E_new=total_energy();
    std::cout<<E_new<<std::endl;
    double L_old=total_L();
    double L_new=total_L();
    int error_times=1;
    int direc=1;
    while (fabs(L_new)>(1e-6)){
        if (((L_new>L_old)and(L_old>0)) or ((L_new>L_old)and(L_old>0))){
            direc*=-1;
            sheer+=direc*s_prim/error_times;
            build_sheer(direc*s_prim/error_times);
        }
        else if (L_new*L_old<0){
            sheer-=direc*s_prim/error_times;
            build_sheer(-1*direc*s_prim/error_times);
            error_times*=2;
        }
        else{
            L_old=L_new;
        }
        sheer+=direc*s_prim/error_times;
        build_sheer(direc*s_prim/error_times);
        L_new=total_L();
        std::cout<<L_new<<" "<<error_times<<" "<<sheer<<std::endl;
    }
    E_new=total_energy();
    std::cout<<E_new<<std::endl;
    std::vector<double> E_L={E_new,L_new};
    return E_L;
}

//读档
double read_world(){
    std::ifstream inFile(filename+".txt");
    std::string line;
    bool need_phys_param=true;
    double L=0;
    count=0;
    std::vector<std::vector<std::vector<Particle>>> savingss(x_num, std::vector<std::vector<Particle>>(y_num, std::vector<Particle>()));
    while (std::getline(inFile, line)) {
            if (need_phys_param){
                std::istringstream iss(line);
                std::string word;
                std::vector<std::string> words;
                while (iss >> word) {
                    words.push_back(word);
                }
                x_tot=std::stoi(words[0]);
                y_tot=std::stoi(words[1]);
                sheer=std::stod(words[2]);
                L=std::stod(words[3]);
                need_phys_param=false;
                double cell=2*max_r+0.1;
                x_num=x_tot/cell;
                y_num=y_tot/cell;
                x_sep=x_tot/double(x_num);
                y_sep=y_tot/double(y_num);
            }
            else{
                std::istringstream iss(line);
                std::string word;
                std::vector<std::string> words;
                while (iss >> word) {
                    words.push_back(word);
                }
                count+=1;
                double x=std::stod(words[0]);
                double y=std::stod(words[1]);
                double r=std::stod(words[2]);
                int id=std::stoi(words[3]);
                Particle particle_i(r,x,y,id);
                savingss[particle_i.x_cell][particle_i.y_cell].push_back(particle_i);
            }
    }
    savings=savingss;
    return L;
}

//开始计算
int main(){
    std::cout<<filename<<" "<<particle_num<<" "<<std::endl;
    double L_0=0;
    L_0=find_md_balance(1e3);
    double E_md=total_energy();
    prdiagnose=true;
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<savings[i][j].size();++k){
                Particle partii=savings[i][j][k];
                place_data[2*partii.id]=partii.x;
                place_data[2*partii.id+1]=partii.y;
            }
        }
    }
    //开始PR平衡
    active_force_vec();
    force_old=force_new;
    std::cout<<force_new[0]<<std::endl;
    double E=10;
    double E_neww=100;
    for (int i=0;i<9;++i){
        int countt=0;
        std::cout<<amp_control<<std::endl;
        while (E/E_neww>(1+amp_control) or (E/E_neww)<(1-amp_control)) {
            E=E_neww;
            countt+=1;
            calculate_direction();
            E_neww=direction_iteration();
            active_force_vec();
        }
        std::cout<<E_neww<<std::endl;
        std::vector<double> E_L_=sheer_itration();
        std::cout<<E_L_[0]<<" "<<E_L_[1]<<std::endl;
        amp_control/=10;
        write_world(E_L_[1]);
    }
    std::cout<<force_new[0]<<std::endl;
}


