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

//信息控件

std::string filename="1006";
std::string testcode="0";

int x_tot=400;
int y_tot=400;

//定义x和y方向格数，总步数计数
int x_num=0;
int y_num=0;
int count=0;
float max_r=1.4;

double passive_force_times=-5;
double active_force_times=0;
double force_beta=0;

double a_f=pow(10,active_force_times/2);
double p_f=pow(10,passive_force_times/2);

int read_txt() {
    std::ifstream inputFile("data.txt");
    if (!inputFile.is_open()) {
        std::cerr << "无法打开文件 data.txt" << std::endl;
    }
    // 按行读取文件内容
    std::string line;
    while (std::getline(inputFile, line)) {
        // 使用stringstream来解析每一行
        std::istringstream iss(line);
        std::string key, value;

        // 假设每行的格式为 key=value
        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            if (key=="active_force_times"){
                active_force_times = std::stod(value);
            } 
            if (key=="force_beta"){
                force_beta = std::stod(value);
            } 
        }
    }
    // 关闭文件
    inputFile.close();

    a_f=pow(10,active_force_times/2);
    p_f=pow(10,passive_force_times/2);

    return 0;
}

int good=read_txt();



//定义弹性系数，阻尼系数
float k=1;
float eta=1;

//定义初时步，形变，预形变
float t_prim=1e-3;
float amp=1e-2;
double sheer=0;

//定义物理界参数
double x_sep=0;
double y_sep=0;



int step_max=1e5;
int heat_random_seed=std::stoi(testcode);
std::normal_distribution<double> normal_dist(0.0, 1.0);
std::random_device rd;
std::mt19937 gen(heat_random_seed);
std::vector<double> time_list={0.01,0.02,0.03,0.05,0.07,
                            0.1,0.2,0.3,0.5,0.7,
                            1,2,3,5,7,
                            10,15,20,30,50,70,
                            1e2,1.5e2,2e2,3e2,5e2,7e2,
                            1e3,1.5e3,2e3,3e3,5e3,7e3,
                            1e4,1.5e4,2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4,
                            1e5,1.5e5,2e5,3e5,4e5,5e5,6e5,7e5,8e5,9e5,
                            1e6,1.5e6,2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6,
                            1e7,1.5e7,2e7,3e7,4e7,5e7,6e7,7e7,8e7,9e7};

//粒子类
class Particle{
public:
    float r;
    double x;
    double y;
    int id;
    int x_cell;
    int y_cell;
    double vx=0;
    double vy=0;
    Particle() : r(0), x(0), y(0), id(0) {}
    Particle(float radius,double particle_x,double particle_y,int particle_id)
        :r(radius),x(particle_x),y(particle_y),id(particle_id){
    x_cell=x/x_sep;
    y_cell=y/y_sep;
    }
};

//创造世界
std::vector<std::vector<std::vector<Particle>>> generate_world() {
    double cell=2*max_r;
    x_num=x_tot/cell;
    y_num=y_tot/cell;
    x_sep=x_tot/double(x_num);
    y_sep=y_tot/double(y_num);
    std::vector<std::vector<std::vector<Particle>>> savings(x_num, std::vector<std::vector<Particle>>(y_num, std::vector<Particle>()));
    return savings;
}

//实例化世界
std::vector<std::vector<std::vector<Particle>>> savings=generate_world();

std::vector<double> final_active_force(80000, 0.0);

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
    std::vector<double> force{force_x,force_y,L_tot};
    return force;
}

//计算移动
std::vector<double> move_to(double fx,double fy,Particle part_self,float t_step){

    double passive_force_x = normal_dist(gen)*pow(t_step,-0.5)*p_f;
    double passive_force_y = normal_dist(gen)*pow(t_step,-0.5)*p_f;

    double active_force_x = normal_dist(gen)*pow(t_step,-0.5)*a_f;
    double active_force_y = normal_dist(gen)*pow(t_step,-0.5)*a_f;

    final_active_force[part_self.id*2]+=force_beta*(active_force_x-final_active_force[part_self.id*2])*t_step;
    final_active_force[part_self.id*2+1]+=force_beta*(active_force_y-final_active_force[part_self.id*2+1])*t_step;

    double x_move=(fx+passive_force_x+final_active_force[part_self.id*2])/eta*t_step;
    double y_move=(fy+passive_force_y+final_active_force[part_self.id*2+1])/eta*t_step;

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

//步进，计算力矩
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
                part_new.vx=move[2];
                part_new.vy=move[3];
                new_savings[part_new.x_cell][part_new.y_cell].push_back(part_new);
            }
        }
    }
    savings=new_savings;
    return sheer_force;
}

//更新自然剪切
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

//保存粒子位置数据
void write_world(int time_step){
    std::ofstream outFile(filename+"_data.txt");
    outFile<<x_tot<<" "<<y_tot<<" "<<sheer<<" "<<time_step<<std::endl;
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<savings[i][j].size();++k){
                Particle part=savings[i][j][k];
                outFile<<part.x<<" "<<part.y<<" "<<part.r<<" "<<part.id<<std::endl;
            }
        }
    }
    outFile.close();
}

//保存剪切数据
void write_force(int time_step,double L_uni,double total_eta,double posi){
    std::ofstream file(filename+"_"+testcode+"_luni.txt", std::ios::app);
    file<<time_step<<" "<<L_uni<<" "<<posi<<" "<<total_eta<<std::endl;
    file.close();
}

//读档
double read_world(){
    std::ifstream inFile(filename+".txt");
    std::string line;
    bool need_phys_param=true;
    double L=0;
    count=0;
    
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
                double cell=2*max_r;
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
                count+=1;
                double x=std::stod(words[0]);
                double y=std::stod(words[1]);
                double r=std::stod(words[2]);
                int id=std::stoi(words[3]);
                Particle particle_i(r,x,y,id);
                savings[particle_i.x_cell][particle_i.y_cell].push_back(particle_i);
            }
    }
    
    for (int i=0;i<x_num;++i) {
        for (int j=0;j<y_num;++j) {
            for (int k=0;k<savings[i][j].size();++k){
                Particle part_old=savings[i][j][k];
            }}}
    
    return L;
}

//开始计算
int main(){
    double total_eta=0;
    double total_eta_posi=0;
    std::ofstream outFile(filename+"_"+testcode+"_luni.txt");
    outFile<<filename<<" "<<testcode<<std::endl;
    outFile.close();
    std::cout<<" "<<filename<<" "<<active_force_times<<" "<<force_beta<<std::endl;
    double L_0=0;
    L_0=read_world();
    sheer+=amp*y_tot;
    build_sheer(amp*y_tot);
    count=0;
    double t_count=0;
    double L_1=step_forward(t_prim);
    L_1=step_forward(t_prim);
    int check_pointer=0;
    double reach_detector=0;
    write_force(0,L_1,0,0);
    for (int i=1;i<200;++i){

        if ((t_count > 1) and (t_prim < 1e-2)){
            t_prim *=1.5849;
        }

        count=0;
        while (count<step_max){
            count+=1;
            t_count+=i*t_prim;
            double L=step_forward(i*t_prim);
            double L_uni=(L-L_0);
            if (L_uni>0){
                total_eta_posi+=L_uni*i*t_prim;
            }
            total_eta+=L_uni*i*t_prim;
            if (t_count>=time_list[check_pointer]){
                std::cout<<time_list[check_pointer]<<" "<<L_uni/L_1<<" "<<total_eta<<" "<<total_eta_posi<<std::endl;
                check_pointer+=1;
                write_force(t_count,L_uni/L_1,total_eta,total_eta_posi);
                write_world(t_count);
            }
        }
    }
}


