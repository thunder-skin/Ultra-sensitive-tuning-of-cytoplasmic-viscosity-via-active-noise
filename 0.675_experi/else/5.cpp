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

//信息控件
double force_times=-3.6;

std::string filename="7005";
std::string testcode="0";

//定义xyz方向格数，总步数计数
int x_num=0;
int y_num=0;
int z_num=0;
int count=0;
float max_r=1.4;

double heat_tense=pow(10,force_times*0.5);

int particle_num=0;

//定义弹性系数，阻尼系数
float k=1;
float eta=1;

//定义初时步，形变，预形变
float t_prim=1e-3;
float amp=1e-2;
double shear=0;

//定义物理界参数
double x_sep=0;
double y_sep=0;
double z_sep=0;

int x_tot=80;
int y_tot=80;
int z_tot=60;

//定义物理世界
std::vector<float> radius={1,1.4};

int step_max=1e3;
int heat_random_seed=std::stoi(testcode);
std::normal_distribution<double> normal_dist(0.0, 1.0);
std::random_device rd;
std::mt19937 gen(heat_random_seed);


//粒子类
class Particle{
public:
    float r;
    double x;
    double y;
    double z;

    int id;
    int x_cell;
    int y_cell;
    int z_cell;

    Particle():r(0),x(0),y(0),z(0),id(0){}
    Particle(float radius,double particle_x,double particle_y,double particle_z,int particle_id)
        :r(radius),x(particle_x),y(particle_y),z(particle_z),id(particle_id){

    x_cell=x/x_sep;
    y_cell=y/y_sep;
    z_cell=z/z_sep;
    }
};

//创造空世界
std::vector<std::vector<std::vector<std::vector<Particle>>>> generate_world() {
    double cell=2*max_r+0.01;
    x_num=x_tot/cell;
    y_num=y_tot/cell;
    z_num=z_tot/cell;
    x_sep=x_tot/double(x_num);
    y_sep=y_tot/double(y_num);
    z_sep=z_tot/double(z_num);
    std::vector<std::vector<std::vector<std::vector<Particle>>>> savings(x_num, std::vector<std::vector<std::vector<Particle>>>(y_num, std::vector<std::vector<Particle>>(z_num,std::vector<Particle>())));
    return savings;
}

//实例化空世界
std::vector<std::vector<std::vector<std::vector<Particle>>>> savings=generate_world();

//计算一对粒子受力
std::vector<double> calculate_force(double x_self,double y_self,double z_self,double x,double y,double z,float r_self,float r,int id1,int id2){
    if (id1==id2){
        std::vector<double> force(4, 0);
        return force;
    }
    else {
        if (y-y_self>2*y_tot/3){
            y_self+=y_tot;
            x_self+=shear;
        }
        else if (y_self-y>2*y_tot/3){
            y_self-=y_tot;
            x_self-=shear;
        }
        if (x_self-x>2*x_tot/3){
            x_self-=x_tot;
        }
        else if (x-x_self>2*x_tot/3){
            x_self+=x_tot;
        }
        if (z_self-z>2*z_tot/3){
            z_self-=z_tot;
        }
        else if (z-z_self>2*z_tot/3){
            z_self+=z_tot;
        }

        if (pow((x_self-x),2)+pow((y_self-y),2)+pow((z_self-z),2)>pow((r_self+r),2)){
            std::vector<double> force(4, 0);
            return force;
        }
        else {
            double delta_r=pow(pow((x_self-x),2)+pow((y_self-y),2)+pow((z_self-z),2),0.5);
            double f_tot=k*(r_self+r-delta_r);
            double fx=f_tot*(x_self-x)/delta_r;
            double fy=f_tot*(y_self-y)/delta_r;
            double fz=f_tot*(z_self-z)/delta_r;
            double Lxy=fx*(y-y_self);
            std::vector<double> force{fx,fy,fz,Lxy};
            return force;
        }
    }
}

//临近受力
std::vector<double> neighbor_force(Particle p){
    double force_x=0;
    double force_y=0;
    double force_z=0;
    double Lxy=0;

    if (p.y_cell==0){
        for (int i=-2;i<=2;++i){
            for (int j=-1;j<=-1;++j){
                for (int k=-1;k<=1;++k){
                    for (int l=0;l<savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)][l];
                        std::vector<double> F_singl=calculate_force(p.x,p.y,p.z,p_e.x,p_e.y,p_e.z,p.r,p_e.r,p.id,p_e.id);
                        force_x+=F_singl[0];
                        force_y+=F_singl[1];
                        force_z+=F_singl[2];
                        Lxy+=F_singl[3];
                    }
                }
            }
        }
        for (int i=-1;i<=1;++i){
            for (int j=0;j<=1;++j){
                for (int k=-1;k<=1;++k){
                    for (int l=0;l<savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)][l];
                        std::vector<double> F_singl=calculate_force(p.x,p.y,p.z,p_e.x,p_e.y,p_e.z,p.r,p_e.r,p.id,p_e.id);
                        force_x+=F_singl[0];
                        force_y+=F_singl[1];
                        force_z+=F_singl[2];
                        Lxy+=F_singl[3];
                    }
                }
            }
        }
    }

    else if (p.y_cell==y_num-1){
        for (int i=-1;i<=1;++i){
            for (int j=-1;j<=0;++j){
                for (int k=-1;k<=1;++k){
                    for (int l=0;l<savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)][l];
                        std::vector<double> F_singl=calculate_force(p.x,p.y,p.z,p_e.x,p_e.y,p_e.z,p.r,p_e.r,p.id,p_e.id);
                        force_x+=F_singl[0];
                        force_y+=F_singl[1];
                        force_z+=F_singl[2];
                        Lxy+=F_singl[3];
                    }
                }
            }
        }
        for (int i=-2;i<=2;++i){
            for (int j=1;j<=1;++j){
                for (int k=-1;k<=1;++k){
                    for (int l=0;l<savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)][l];
                        std::vector<double> F_singl=calculate_force(p.x,p.y,p.z,p_e.x,p_e.y,p_e.z,p.r,p_e.r,p.id,p_e.id);
                        force_x+=F_singl[0];
                        force_y+=F_singl[1];
                        force_z+=F_singl[2];
                        Lxy+=F_singl[3];
                    }
                }
            }
        }
    }

    else{ 
        for (int i=-1;i<=1;++i){
            for (int j=-1;j<=1;++j){
                for (int k=-1;k<=1;++k){
                    for (int l=0;l<savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)][l];
                        std::vector<double> F_singl=calculate_force(p.x,p.y,p.z,p_e.x,p_e.y,p_e.z,p.r,p_e.r,p.id,p_e.id);
                        force_x+=F_singl[0];
                        force_y+=F_singl[1];
                        force_z+=F_singl[2];
                        Lxy+=F_singl[3];
                    }
                }
            }
        }
    }

    std::vector<double> force{force_x,force_y,force_z,Lxy};
    return force;
}

//计算一对粒子能量
double calculate_energy(double x_self,double y_self,double z_self,double x,double y,double z,float r_self,float r,int id1,int id2){
    if (id1==id2){
        return 0;
    }
    else {
        if (y-y_self>2*y_tot/3){
            y_self+=y_tot;
            x_self+=shear;
        }
        else if (y_self-y>2*y_tot/3){
            y_self-=y_tot;
            x_self-=shear;
        }
        if (x_self-x>2*x_tot/3){
            x_self-=x_tot;
        }
        else if (x-x_self>2*x_tot/3){
            x_self+=x_tot;
        }
        if (z_self-z>2*z_tot/3){
            z_self-=z_tot;
        }
        else if (z-z_self>2*z_tot/3){
            z_self+=z_tot;
        }
        if (pow((x_self-x),2)+pow((y_self-y),2)+pow((z_self-z),2)>pow((r_self+r),2)) {
            return 0;
        }
        else {
            double delta_r=pow(pow((x_self-x),2)+pow((y_self-y),2)+pow((z_self-z),2),0.5);
            double E=0.5*pow((r+r_self-delta_r),2);
            return E;
        }
    }
}

//临近能量
double neighbor_energy(Particle p){
    double E_tot=0;

    if (p.y_cell==0){
        for (int i=-2;i<=2;++i){
            for (int j=-1;j<=-1;++j){
                for (int k=-1;k<=1;++k){
                    for (int l=0;l<savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)][l];
                        double E_singl=calculate_energy(p.x,p.y,p.z,p_e.x,p_e.y,p_e.z,p.r,p_e.r,p.id,p_e.id);
                        E_tot+=E_singl;
                    }
                }
            }
        }
        for (int i=-1;i<=1;++i){
            for (int j=0;j<=1;++j){
                for (int k=-1;k<=1;++k){
                    for (int l=0;l<savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)][l];
                        double E_singl=calculate_energy(p.x,p.y,p.z,p_e.x,p_e.y,p_e.z,p.r,p_e.r,p.id,p_e.id);
                        E_tot+=E_singl;
                    }
                }
            }
        }
    }

    else if (p.y_cell==y_num-1){
        for (int i=-1;i<=1;++i){
            for (int j=-1;j<=0;++j){
                for (int k=-1;k<=1;++k){
                    for (int l=0;l<savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)][l];
                        double E_singl=calculate_energy(p.x,p.y,p.z,p_e.x,p_e.y,p_e.z,p.r,p_e.r,p.id,p_e.id);
                        E_tot+=E_singl;
                    }
                }
            }
        }
        for (int i=-2;i<=2;++i){
            for (int j=1;j<=1;++j){
                for (int k=-1;k<=1;++k){
                    for (int l=0;l<savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)][l];
                        double E_singl=calculate_energy(p.x,p.y,p.z,p_e.x,p_e.y,p_e.z,p.r,p_e.r,p.id,p_e.id);
                        E_tot+=E_singl;
                    }
                }
            }
        }
    }

    else{ 
        for (int i=-1;i<=1;++i){
            for (int j=-1;j<=1;++j){
                for (int k=-1;k<=1;++k){
                    for (int l=0;l<savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[fmod(p.x_cell+i+x_num,x_num)][fmod(p.y_cell+j+y_num,y_num)][fmod(p.z_cell+k+z_num,z_num)][l];
                        double E_singl=calculate_energy(p.x,p.y,p.z,p_e.x,p_e.y,p_e.z,p.r,p_e.r,p.id,p_e.id);
                        E_tot+=E_singl;
                    }
                }
            }
        }
    }

    return E_tot;
}

//计算总能量
double total_energy(){
    double Energy_total=0;
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<z_num;++k) {
                for (int l=0;l<savings[i][j][k].size();++l){
                    Energy_total+=neighbor_energy(savings[i][j][k][l]);
                }
            }
        }
    }
    return Energy_total;
}

//计算MD移动
std::vector<double> move_to(double fx,double fy,double fz,Particle part_self,float t_step){

    double ra_x = normal_dist(gen)*pow(t_step,-0.5)*heat_tense;
    double ra_y = normal_dist(gen)*pow(t_step,-0.5)*heat_tense;
    double ra_z = normal_dist(gen)*pow(t_step,-0.5)*heat_tense;
    
    double x_move=(fx+ra_x)/eta*t_step;
    double y_move=(fy+ra_y)/eta*t_step;
    double z_move=(fz+ra_z)/eta*t_step;

    double x_alter=0;
    if (y_move+part_self.y>y_tot){
        x_alter=-shear;
    }
    else if (y_move+part_self.y<0){
        x_alter=shear;
    }

    std::vector<double> move{fmod(part_self.x+x_move+x_alter+x_tot,x_tot),fmod(part_self.y+y_move+y_tot,y_tot),fmod(part_self.z+z_move+z_tot,z_tot)};
    return move;
}

//步进计算力矩
double step_forw(float t_step){
    std::vector<std::vector<std::vector<std::vector<Particle>>>> new_savings(x_num, std::vector<std::vector<std::vector<Particle>>>(y_num, std::vector<std::vector<Particle>>(z_num,std::vector<Particle>())));
    double shear_force=0;
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<z_num;++k){
                for (int l=0;l<savings[i][j][k].size();++l){
                    Particle part_old=savings[i][j][k][l];
                    std::vector<double> force=neighbor_force(part_old);
                    std::vector<double> move=move_to(force[0],force[1],force[2],part_old,t_step);
                    shear_force+=force[3];
                    Particle part_new(part_old.r,move[0],move[1],move[2],part_old.id);
                    new_savings[part_new.x_cell][part_new.y_cell][part_new.z_cell].push_back(part_new);
                }
            }
        }
    }
    savings=new_savings;
    return shear_force;
}

//更新剪切
void build_shear(double shear){
    std::vector<std::vector<std::vector<std::vector<Particle>>>> new_savings(x_num, std::vector<std::vector<std::vector<Particle>>>(y_num, std::vector<std::vector<Particle>>(z_num,std::vector<Particle>())));
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<z_num;++k) {
                for (int l=0;l<savings[i][j][k].size();++l){
                    Particle part_old=savings[i][j][k][l];
                    Particle part_new(part_old.r,fmod(part_old.x+part_old.y*shear/y_tot+x_tot,x_tot),part_old.y,part_old.z,part_old.id);
                    new_savings[part_new.x_cell][j][k].push_back(part_new);
                }
            }
        }
    }
    savings=new_savings;
}

//保存粒子位置数据
void write_world(double time,double t_step){
    std::ofstream outFile(filename+"_data.txt");
    outFile<<x_tot<<" "<<y_tot<<" "<<z_tot<<" "<<particle_num<<" "<<shear<<" "<<time<<" "<<t_step<<std::endl;
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<z_num;++k) {
                for (int l=0;l<savings[i][j][k].size();++l){
                    Particle part=savings[i][j][k][l];
                    std::ostringstream oss;
                    oss<<std::fixed<<std::setprecision(14);
                    oss<<part.x<<" "<<part.y<<" "<<part.z<<" "<<part.r<<" "<<part.id<<std::endl;
                    outFile<<oss.str();
                }
            }
        }
    }
    outFile.close();
}

//保存剪切数据
void write_force(double time_step,double L_uni,double total_eta,double posi,double E){
    std::ofstream file(filename+"_"+testcode+"_luni.txt", std::ios::app);
    std::ostringstream oss;
    oss<<std::fixed<<std::setprecision(7);
    oss<<time_step<<" "<<L_uni<<" "<<posi<<" "<<total_eta<<" "<<E<<std::endl;
    file<<oss.str();
    file.close();
}

//读档
double read_world(){
    std::ifstream inFile(filename+".txt");
    std::string line;
    bool need_phys_param=true;
    double L=0;
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
            z_tot=std::stoi(words[2]);
            particle_num=std::stoi(words[3]);
            shear=std::stod(words[4]);
            L=std::stod(words[5]);

            need_phys_param=false;
            double cell=2*max_r+0.01;

            x_num=x_tot/cell;
            y_num=y_tot/cell;
            x_sep=x_tot/double(x_num);
            y_sep=y_tot/double(y_num);

            std::vector<std::vector<std::vector<Particle>>> savings(x_num, std::vector<std::vector<Particle>>(y_num, std::vector<Particle>()));
        }
        else{
            std::istringstream iss(line);
            std::string word;
            std::vector<std::string> words;
            while (iss >> word) {
                words.push_back(word);
            }
            double x=std::stod(words[0]);
            double y=std::stod(words[1]);
            double z=std::stod(words[2]);
            double r=std::stod(words[3]);
            int id=std::stoi(words[4]);
            Particle particle_i(r,x,y,z,id);
            savings[particle_i.x_cell][particle_i.y_cell][particle_i.z_cell].push_back(particle_i);
        }
    }
    return L;
}

//开始计算
int main(){

    //读取世界，施加初始剪切
    double L_0=read_world();
    shear+=amp*y_tot;
    build_shear(amp*y_tot);

    //定义一些需要记录的参数
    double total_eta=0;
    double total_eta_posi=0;
    double t_count=0;

    //初期弛豫，减少涨落
    double L_1=0;
    for (int i=1;i<10;++i){
        L_1=step_forw(t_prim);
    }
    
    //记录初始数据
    std::cout<<filename<<" "<<force_times<<" "<<testcode<<" "<<L_1<<std::endl;
    std::ofstream outFile(filename+"_"+testcode+"_luni.txt");
    outFile<<filename<<" "<<force_times<<" "<<testcode<<" "<<L_1<<std::endl;
    outFile.close();

    std::vector<double> time_list;
    for (double x = -1.5; x <= 7; x += 0.05) {
        time_list.push_back(std::pow(10, x));
    }
    
    int check_pointer=0;
    double reach_detector=0;
    double t_stepp=t_prim;
    double step_count=0;
    double L_uni=0;

    while (t_count<1e7){
        step_count=0;
        if (t_stepp<1e-1){
            t_stepp*=1.05;
        }
        else{
            t_stepp*=1.001;
        }
        while (step_count<step_max){
            step_count+=1;
            t_count+=t_stepp;
            L_uni=step_forw(t_stepp)-L_0;
            if (L_uni>0){
                total_eta_posi+=L_uni*t_stepp;
            }
            total_eta+=L_uni*t_stepp;
            if (t_count>=time_list[check_pointer]){
                double E=total_energy();
                std::cout<<time_list[check_pointer]<<" "<<L_uni/L_1<<" "<<total_eta<<" "<<t_stepp<<" "<<E<<std::endl;
                check_pointer+=1;
                write_force(t_count,L_uni/L_1,total_eta,total_eta_posi,E);
                write_world(t_count,t_stepp);

            }
        }
    }
}


