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
int z_num=0;
int count=0;
std::string filename="1002";
int particle_num=553;
double ccount_max=100;

//定义物理界参数
double x_sep=0;
double y_sep=0;
double z_sep=0;

int x_tot=10;
int y_tot=10;
int z_tot=10;

float k=1;
float eta=1;

double t_prim=1e-2;
double shear=0;

//定义工作
bool prdiagnose=false;
std::vector<double> force_new(3*particle_num,0);
std::vector<double> force_old(3*particle_num,0);
std::vector<double> velocy_new(3*particle_num,0);
std::vector<double> velocy_old(3*particle_num,0);
std::vector<double> place_data(3*particle_num,0);

//定义物理世界
std::vector<float> radius={1};
float max_r=1;
int balance_max=1e4;
int random_seed=std::stoi(filename);
std::mt19937 gen(random_seed);
std::random_device rd;

int m(int a,int b){
    if (a>=b){
        a-=b;
    }
    if (a>=b){
        a-=b;
    }
    if (a<0){
        a+=b;
    }
    if (a<0){
        a+=b;
    }
    return a;
}

double mm(double a,double b){
    while (a>b){
        a-=b;
    } 
    while (a<0){
        a+=b;
    }
    return a;
}

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

//创造随机初始世界
std::vector<std::vector<std::vector<std::vector<Particle>>>> generate_world() {
    double cell=2*max_r+0.01;
    x_num=x_tot/cell;
    y_num=y_tot/cell;
    z_num=z_tot/cell;
    x_sep=x_tot/double(x_num);
    y_sep=y_tot/double(y_num);
    z_sep=z_tot/double(z_num);
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    std::vector<std::vector<std::vector<std::vector<Particle>>>> savings(x_num, std::vector<std::vector<std::vector<Particle>>>(y_num, std::vector<std::vector<Particle>>(z_num,std::vector<Particle>())));
    for (int i=0;i<particle_num;++i) {
        double random_num=0;
        random_num=dis(gen);
        double r_i=0;
        
        r_i=log(1.1052+(2.7182-1.1052)*random_num);
        random_num=dis(gen);
        double x_i=x_tot*random_num;
        random_num=dis(gen);
        double y_i=y_tot*random_num;
        random_num=dis(gen);
        double z_i=z_tot*random_num;
        Particle particle_i(r_i,x_i,y_i,z_i,i);
        savings[particle_i.x_cell][particle_i.y_cell][particle_i.z_cell].push_back(particle_i);
    }
    return savings;
}

//实例化随机初始世界
std::vector<std::vector<std::vector<std::vector<Particle>>>> savings=generate_world();

//计算一对粒子受力
std::vector<double> calculate_force(double x_self,double y_self,double z_self,double x,double y,double z,float r_self,float r,int id1,int id2){
    if (id1==id2){
        std::vector<double> force(3, 0);
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
            std::vector<double> force(3, 0);
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
                    for (int l=0;l<savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)][l];
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
                    for (int l=0;l<savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)][l];
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
                    for (int l=0;l<savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)][l];
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
                    for (int l=0;l<savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)][l];
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
                    for (int l=0;l<savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)][l];
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
    
    if (prdiagnose==true){
        force_new[3*p.id]=force_x;
        force_new[3*p.id+1]=force_y;
        force_new[3*p.id+2]=force_z;
    }

    std::vector<double> force{force_x,force_y,force_z,Lxy};
    return force;
}

//临近能量
double neighbor_energy(Particle p){
    double E_tot=0;

    if (p.y_cell==0){
        for (int i=-2;i<=2;++i){
            for (int j=-1;j<=-1;++j){
                for (int k=-1;k<=1;++k){
                    for (int l=0;l<savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)][l];
                        double E_singl=calculate_energy(p.x,p.y,p.z,p_e.x,p_e.y,p_e.z,p.r,p_e.r,p.id,p_e.id);
                        E_tot+=E_singl;
                    }
                }
            }
        }
        for (int i=-1;i<=1;++i){
            for (int j=0;j<=1;++j){
                for (int k=-1;k<=1;++k){
                    for (int l=0;l<savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)][l];
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
                    for (int l=0;l<savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)][l];
                        double E_singl=calculate_energy(p.x,p.y,p.z,p_e.x,p_e.y,p_e.z,p.r,p_e.r,p.id,p_e.id);
                        E_tot+=E_singl;
                    }
                }
            }
        }
        for (int i=-2;i<=2;++i){
            for (int j=1;j<=1;++j){
                for (int k=-1;k<=1;++k){
                    for (int l=0;l<savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)][l];
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
                    for (int l=0;l<savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)].size();++l){
                        Particle p_e=savings[m(p.x_cell+i+x_num,x_num)][m(p.y_cell+j+y_num,y_num)][m(p.z_cell+k+z_num,z_num)][l];
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

//计算总剪切矩
double total_L(){
    double L=0;
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<z_num;++k) {
                for (int l=0;l<savings[i][j][k].size();++l){
                    Particle part_old=savings[i][j][k][l];
                    std::vector<double> force=neighbor_force(part_old);
                    L+=force[3];
                }
            }
        }
    }
    return L;
}

//存档
void write_world(){
    std::string name=std::to_string(random_seed);
    std::ofstream outFile(filename+".txt");
    double L=total_L();
    outFile<<x_tot<<" "<<y_tot<<" "<<z_tot<<" "<<particle_num<<" "<<shear<<" "<<L<<std::endl;
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<z_num;++k) {
                for (int l=0;l<savings[i][j][k].size();++l){
                    Particle part=savings[i][j][k][l];
                    std::ostringstream oss;
                    oss<<std::fixed<<std::setprecision(17);
                    oss<<part.x<<" "<<part.y<<" "<<part.z<<" "<<part.r<<" "<<part.id<<std::endl;
                    outFile<<oss.str();
                }
            }
        }
    }
    outFile.close();
}

//计算MD移动
std::vector<double> move_to(double fx,double fy,double fz,Particle part_self,float t_step){
    double x_move=fx/eta*t_step;
    double y_move=fy/eta*t_step;
    double z_move=fz/eta*t_step;
    double x_alter=0;
    if (y_move+part_self.y>y_tot){
        x_alter=-shear;
    }
    else if (y_move+part_self.y<0){
        x_alter=shear;
    }
    std::vector<double> move{mm(part_self.x+x_move+x_alter,x_tot),mm(part_self.y+y_move,y_tot),mm(part_self.z+z_move,z_tot)};
    return move;
}

//MD步进
void step_forward(float t_step){
    std::vector<std::vector<std::vector<std::vector<Particle>>>> new_savings(x_num, std::vector<std::vector<std::vector<Particle>>>(y_num, std::vector<std::vector<Particle>>(z_num,std::vector<Particle>())));
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<z_num;++k){
                for (int l=0;l<savings[i][j][k].size();++l){
                    Particle part_old=savings[i][j][k][l];
                    std::vector<double> force=neighbor_force(part_old);
                    std::vector<double> move=move_to(force[0],force[1],force[2],part_old,t_step);
                    Particle part_new(part_old.r,move[0],move[1],move[2],part_old.id);
                    new_savings[part_new.x_cell][part_new.y_cell][part_new.z_cell].push_back(part_new);
                }
            }
        }
    }
    savings=new_savings;
    new_savings.clear(); // 清除所有元素
    std::vector<std::vector<std::vector<std::vector<Particle>>>>().swap(new_savings);
}

//建立MD平衡
void find_md_balance(int md_max){
    int count=0;
    while (count<=md_max){
        count+=1;
        step_forward(1e-2);
    }
    double E_2=total_energy();
}

//开始计算
int main(){
    std::cout<<filename<<" "<<particle_num<<" "<<std::endl;

    write_world();
    double L_0=0;
    find_md_balance(1000);
    double ccount=0;
    double E_0=total_energy();
    double E_new=E_0;
    while (ccount<ccount_max){
        std::cout<<ccount<<" "<<E_new<<std::endl;
        ccount+=1;
        find_md_balance(1000);
        E_new=total_energy();
        write_world();
    }
}


