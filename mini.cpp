#include <iostream>
#include <string>
#include <bits/stdc++.h>
using namespace std;
#define maxsize 1024
#define epsilon 0.707
#define image_pad 5
#define corner 10
#define c_1 1 
#define smoothness 2

// global volume of voxels
std::map<int,std::map<int,std::bitset<2*maxsize> > > boolarr;


typedef struct{
	int x;
	int y;
	int z;
	int index;
	bool keep;
}points;


typedef struct normals{
	double n_x;
	double n_y;
	double n_z;

	normals(double a,double b,double c):n_x(a),n_y(b),n_z(c){}
	normals():n_x(0),n_y(0),n_z(0){}
}normals;
// point structure for queue
typedef struct point{
	int x;
	int y;
	int z;
	
	point(int a,int b,int c):x(a),y(b),z(c){}
	point():x(0),y(0),z(0){}
}point;

double dist(double a,double b, double c, double d)
{

	return sqrt((a-b)*(a-b) + (c-d)*(c-d));
}


double dot_prod(normals first, normals last){
	return (first.n_x*last.n_x + first.n_y*last.n_y + first.n_z*last.n_z)/((sqrt(first.n_x*first.n_x + first.n_y*first.n_y + first.n_z*first.n_z))*(sqrt(last.n_x*last.n_x + last.n_y*last.n_y + last.n_z*last.n_z)));
}

struct by_keep { 
    bool operator()(points const &a, points const &b) { 
        return a.keep>b.keep;
    }
};


struct by_x { 
    bool operator()(points const &a, points const &b) { 
        return a.x<b.x;
    }
};

// takes too much memory
/*
void floodfill(int x, int y, int z){
	// corner cases
	if(x>=maxsize || x<0 || y>=maxsize || y<0 || z>=maxsize || z<0)
		return;

	cout<<x<<" "<<y<<" "<<z;
	
	if(boolarr[x][y].test(2*z) && boolarr[x][y].test(2*z+1)){
		cout<<" 11\n";
	}
	else if(boolarr[x][y].test(2*z) && !boolarr[x][y].test(2*z+1)){
		
		cout<<" 10\n";
	}
	else if(!boolarr[x][y].test(2*z) && boolarr[x][y].test(2*z+1))
		cout<<" 01\n";
	else
		cout<<" 00\n";
	// surface --> 11, internal ----> 10
	if(boolarr[x][y].test(2*z))
	{
		boolarr[x][y].set(2*z+1);
		return;
	}
	// external ---> 01
	else if(boolarr[x][y].test(2*z+1))
	{
		return;
	}

	// external point ----> 01
	boolarr[x][y].set(2*z+1);	

	// 6-neighbourhood movement
	
	floodfill(x+1,y,z);
	
	cout<<"top\n";
	floodfill(x,y+1,z);
	
	cout<<"front\n";
	floodfill(x,y,z+1);
	cout<<"left\n";
	floodfill(x-1,y,z);
	cout<<"bottom\n";
	floodfill(x,y-1,z);
	cout<<"back\n";
	floodfill(x,y,z-1);
	return;
}
*/



void floodfill(point start,int maxx,int maxy,int maxz){
 	std::stack<point> stack;
 	int minx = start.x;
 	int miny = start.y;
 	int minz = start.z;

    stack.push(start);
    while (!stack.empty())
    {
        point p = stack.top();
        stack.pop();
        int x = p.x;
        int y = p.y;
        int z = p.z;
        //cout<<x<<" "<<y<<" "<<z<<"\n";
        if (y < miny || y > maxy - 1 || x < minx || x > maxx - 1 || z<minz || z>maxz - 1)
            continue;
        if(boolarr[x][y].test(2*z))
		{
			boolarr[x][y].set((2*z)+1);
			continue;
		}
		// external ---> 01
		else if(boolarr[x][y].test(2*z+1))
		{
			continue;
		}
		boolarr[x][y].set((2*z)+1);	

        if(!(x+1 > maxx - 1 ) && !(!boolarr[x+1][y].test(2*z) && boolarr[x+1][y].test(2*z+1)))
        	stack.push(point(x+1,y,z));
        if(!(y+1 > maxy - 1 ) && !(!boolarr[x][y+1].test(2*z) && boolarr[x][y+1].test(2*z+1)) )
        	stack.push(point(x,y+1,z));
        if(!(z+1 > maxz - 1 ) && !(!boolarr[x][y].test(2*(z+1)) && boolarr[x][y].test(2*(z+1)+1)) )
        	stack.push(point(x,y,z+1));
        if(!(x-1 < minx ) && !(!boolarr[x-1][y].test(2*z) && boolarr[x-1][y].test(2*z+1)) )
        	stack.push(point(x-1,y,z));
        if(!(y-1 < miny ) && !(!boolarr[x][y-1].test(2*z) && boolarr[x][y-1].test(2*z+1)) )
        	stack.push(point(x,y-1,z));
        if(!(z-1 < minz ) && !(!boolarr[x][y].test(2*(z-1)) && boolarr[x][y].test(2*(z-1)+1)))
       		stack.push(point(x,y,z-1));

    }
}

void hsv2rgb(double h,double s,double v,int *R,int *G,int *B)
{

	if(h>360 || h<0.0)
	{
		*R = (int)255*(v-v*s);
		*G = (int)255*(v-v*s);
		*B = (int)255*(v-v*s);

	}

	double c = v*s;
	double R1=0.0,G1=0.0,B1=0.0;
	double Hprime = h/60;
	double m = v-c;
	int temp = ((int)Hprime)%2 - 1;
	if(temp<0)
		temp = -temp;

	double x = c * (1 - temp);
	if(Hprime <= 1.0)
	{
		R1 = c;
		G1 = x;
	}
	else if(Hprime<=2.0)
	{
		R1 = x;
		G1 = c;
	}

	else if(Hprime <= 3.0)
	{
		G1 = c;
		B1 = x;
	}
	else if(Hprime<=4.0)
	{
		G1 = x;
		B1 = c;
	}

	else if(Hprime <= 5.0)
	{
		B1 = c;
		R1 = x;
	}
	else if(Hprime<=6.0)
	{
		B1 = x;
		R1 = c;
	}

	*R = (int)255*(R1+m);
	*G = (int)255*(G1+m);
	*B = (int)255*(B1+m);
}


void final_step(double **present,double mininten,double maxinten,int x_left, int x_right, int y_left, int y_right){
	double Rprime,Gprime,Bprime,Cmax,Cmin,del,S,V,H;
	long long x,R,G,B;
	int checker,i,j;
	while(1)
	{
		cout<<"Enter a color in hexadecimal format\n";
		cin>> hex >> x;
		//cout<< x;
		if(x>0xffffff){
			cout<<"Invalid Input.\n";
		}
		else
			break;
	}

	B = x%256;
	x = x/256;
	G = x%256;
	x = x/256;
	R = x%256;

	Rprime = R/255.0;
	Gprime = G/255.0;
	Bprime = B/255.0;

	
	if(Rprime >= Gprime and Rprime>=Bprime)
	{
		Cmax = Rprime;
		checker = 0;
	}
	else if(Gprime >= Bprime)
	{
		Cmax = Gprime;
		checker = 1;
	}
	else 
	{
		Cmax = Bprime;
		checker = 2;
	}	


	if(Rprime <= Gprime and Rprime <= Bprime)
		Cmin = Rprime;
	else if(Gprime <= Bprime)
		Cmin = Gprime;
	else Cmin = Bprime;

	del = Cmax-Cmin;

	V =Cmax;
	if(Cmax != 0.0)
		S = del/Cmax;
	else
		S = 0.0;

	if(del != 0.0)
	{
		if(checker==0)
			H = ((Gprime-Bprime)/del);
		else if(checker == 1)
			H = 2 + (Bprime-Rprime)/del;
		else
			H = 4 + (Rprime-Gprime)/del;
		H *= 60;
		if(H<0)
			H+=360.0;

	}
	else
		H = 1000.0;
	

	ofstream myfile;
	myfile.open ("output1.svg");

	myfile<<"<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n";
	myfile<<"<rect x=\"0\" y=\"0\" width=\""<<(x_right-x_left+1)<<"\" height=\""<<(y_right-y_left+1)<<"\" stroke=\"black\" stroke-width=\"0.3\" fill=\"#000000\" />\n";

	for(i=x_left;i<=x_right;i++)
		for(j=y_left;j<=y_right;j++)
		{
			if(present[i-x_left][j-y_left]==0.0)
			{
				myfile<<"<rect x=\""<< (i-x_left) <<"\" y=\""<<  (j-y_left) <<"\" width=\"8\" height=\"8\" style=\"fill:#ffffff;fill-opacity:0.5\" />\n";
				continue;
			}


			int R_new,G_new,B_new;
			double v_val;
			v_val = (present[i-x_left][j-y_left] - mininten)/(maxinten - mininten);
			hsv2rgb(H,S,v_val,&R_new,&G_new,&B_new);

			//cout << H <<"	"<<S<<"	"<<v_val<<"\n";
			//cout << R_new <<"	"<< G_new <<"	"<< B_new<<"\n";
			std::stringstream stream;
			stream << setfill('0') << setw(2) << std::hex << R_new;
			string R_res(stream.str());

			stream  << setfill('0') << setw(2) << std::hex << G_new;
			string G_res(stream.str());

			stream  << setfill('0') << setw(2) << std::hex << B_new;
			string B_res(stream.str());

			//cout<<B_res<<"\n";

			myfile<<"<rect x=\""<< (i-x_left) <<"\" y=\""<<  (j-y_left) <<"\" width=\"8\" height=\"8\" style=\"fill:#"<<B_res<<";fill-opacity:0.5\" />\n";
		}

	myfile<<"</svg>\n"; 
	myfile.close();	

}


int  main()
{
	// variables
	int n, a, b, c, i, j, k, l, minx = maxsize, miny = maxsize, minz = maxsize,tot = 0, maxx = -maxsize, maxy = -maxsize, maxz = -maxsize, temp_z;
	string line;
	double minframe_x,minframe_y,maxframe_x,maxframe_y,max_intens,min_intens;
	int x_left,y_left,x_right,y_right;
	double **present;
	double K_a,I_a,I_p,K_s,K_d;
	//std::map<int,std::map<int,std::bitset<2*maxsize> > > boolarr;

	int x_v,y_v,z_v,x_c,y_c,z_c,x_i,y_i,z_i;
	double temp_x,temp_y;
	std::map<int,std::map<int,std::map<int,normals> > > check;
	// array of points
	std::vector<points> pointarr;	
	
	ifstream file1;
	file1.open("circle.txt");

	//for(i=0;i<maxsize;i++)
	//	for(j=0;j<maxsize;j++)
	//		boolarr[i][j].reset(); 	

	//cout<<"Now enter the coordinates\n";
	// read each and every point

	if(file1.is_open())
	{
		while(getline(file1,line))
		{	
			//cout<<line<<"\n";
			if(tot==0)
				{
					n = atoi(line.c_str());
					tot+=1;
				}
			else
			{
				std::stringstream stream(line);
				//centres temp;

				for(i=0;i<3;i++){
					
					if(i==0)
						stream >> a;
					else if(i==1)
						stream >> b;
					else
						stream >> c;
				}

				points point;
	
				point.x = a;
				point.y = b;
				point.z = c;

				point.index = i;
				pointarr.push_back(point);

				if(minx > a)
					minx = a;
				if(miny > b)
					miny = b;
				if(minz > c)
					minz = c;
				if(maxx < a)
					maxx = a;
				if(maxy < b)
					maxy = b;
				if(maxz < c)
					maxz = c;

				
			}
		}
	}
	file1.close();
	



	cout<<minx<<" "<<miny<<" "<<minz<<"\n";
	cout<<maxx<<" "<<maxy<<" "<<maxz<<"\n";


	// make all non negative 
	for(i=0;i<n;i++)
	{	
		pointarr[i].x = pointarr[i].x - minx + corner;
		pointarr[i].y = pointarr[i].y - miny + corner;
		pointarr[i].z = pointarr[i].z - minz + corner;

		// maintain the 3D array of bits
		// cout<<pointarr[i].x<<" "<<pointarr[i].y<<" "<<pointarr[i].z<<"\n";
		(boolarr[pointarr[i].x][pointarr[i].y]).set(2*(pointarr[i].z)); 
		//check[pointarr[i].x][pointarr[i].y][pointarr[i].z] = false;
	}

	// cout<<boolarr[82][45].to_string<char,std::string::traits_type,std::string::allocator_type>();
	

	
	maxz -= (minz - corner);
	maxy -= (miny - corner);
	maxx -= (minx - corner);
	
	cout<<"Enter view position\n";
	cin >> x_v >> y_v >> z_v;

	cout<<"Enter camera position\n";
	cin >> x_c >> y_c >> z_c;


	//cout<<"Enter left hand corner of the image plane\n";
	//cin >> x_i >> y_i >> z_i;

	
	// tranform the vertices since we have translated all the voxels
	x_v -= (minx-corner);
	y_v -= (miny-corner);
	z_v -= (minz-corner);

	x_c -= (minx-corner);
	y_c -= (miny-corner);
	z_c -= (minz-corner);

	

	minx = corner;
	miny = corner;
	minz = corner;

	// translate the entire object space to z = -corner
	
	z_i = maxz + corner;


	// find the normals at each point
	for(i=0;i<n;i++)
	{

		int countx = 0, county = 0, countz = 0;
		// (boolarr[pointarr[i].x][pointarr.y]).test(pointarr[i].z);
		// x gradient (projection on yz) 


		for(k = std::max<int>(pointarr[i].y-1,0); k <= std::min<int>(pointarr[i].y+1,maxsize); k++) 
			for(j = std::max<int>(pointarr[i].z-1,0); j <= std::min<int>(pointarr[i].z+1,maxsize); j++)
			{
				a = std::max<int>(pointarr[i].x-1,0);
				b = a + 1;
				c = std::min<int>(b + 1, maxsize);

				if(boolarr[a][k].test(2*j) || (boolarr[b][k].test(2*j) && k!=pointarr[i].y && j!=pointarr[i].z) || boolarr[c][k].test(2*j))
					countx++;
			}



		// y gradient (projection on xz) 
		for(k = std::max<int>(pointarr[i].x-1,0); k <= std::min<int>(pointarr[i].x+1,maxsize); k++) 
			for(j = std::max<int>(pointarr[i].z-1,0); j <= std::min<int>(pointarr[i].z+1,maxsize); j++)
			{
				a = std::max<int>(pointarr[i].y-1,0);
				b = a + 1;
				c = std::min<int>(b + 1, maxsize);

				if(boolarr[k][a].test(2*j) || (boolarr[k][b].test(2*j) && k!=pointarr[i].x && j!=pointarr[i].z) || boolarr[k][c].test(2*j))
					county++;
			}	


		// z  gradeint (projection on xz) 
		for(k = std::max<int>(pointarr[i].x-1,0); k <= std::min<int>(pointarr[i].x+1,maxsize); k++) 
			for(j = std::max<int>(pointarr[i].y-1,0); j <= std::min<int>(pointarr[i].y+1,maxsize); j++)
			{
				a = std::max<int>(pointarr[i].z-1,0);
				b = a + 1;
				c = std::min<int>(b + 1, maxsize);

				if(boolarr[k][j].test(2*a) || (boolarr[k][j].test(2*b) && k!=pointarr[i].x && j!=pointarr[i].y) || boolarr[k][j].test(2*c))
					countz++;
			}	

		if(countx == 0 && county==0 && countz == 0)
		{
			countx = 1;
			county = 1;
			countz = 1;
		}	

		double sol = sqrt(countx*countx + county*county + countz*countz);
		//pointarr[i].n_x = countx/sol;
		//pointarr[i].n_y = county/sol;
		//pointarr[i].n_z = countz/sol;
		check[pointarr[i].x][pointarr[i].y][pointarr[i].z] = normals(countx/sol,county/sol,countz/sol);
		
	}

	

	
	//for(i=0;i<n;i++){
	//	cout<<pointarr[i].x<<" "<<pointarr[i].y<<" "<<pointarr[i].z<<"\n";
	//}


	sort(pointarr.begin(),pointarr.end(),by_x());

	// cout<<maxx<<" "<<maxy<<" "<<maxz<<"\n";

	//cout<<"\n"<<boolarr[82][45].to_string<char,std::string::traits_type,std::string::allocator_type>();


	// flood fill
	floodfill(point(minx-1,miny-1,minz-1),min(maxx+2,maxsize),min(maxy+2,maxsize),min(maxz+2,maxsize));

	// cout<<"\n"<<boolarr[82][45].to_string<char,std::string::traits_type,std::string::allocator_type>();

	
	
	// sign checking of normals
	for(i=0;i<n;i++){

		int checkx = 0,checky=0,checkz=0;
		int tx,ty,tz;
		// check for normal x
		a = pointarr[i].x;
		b = pointarr[i].y;
		c = pointarr[i].z;

		// check if (a-1,b,c) is internal(*0) and (a+1,b,c) is external (01) or viceversa
		if(!boolarr[a-1][b].test(2*c+1) && (!boolarr[a+1][b].test(2*c) && boolarr[a+1][b].test(2*c+1)))
		{
			checkx = 1;
		}
		else if(!boolarr[a+1][b].test(2*c+1) && (!boolarr[a-1][b].test(2*c) && boolarr[a-1][b].test(2*c+1)))
		{
			checkx = -1;	
		}

		// check if (a,b-1,c) is internal(*0) and (a,b+1,c) is external(01) or viceversa
		if(!boolarr[a][b-1].test(2*c+1) && (!boolarr[a][b+1].test(2*c) && boolarr[a][b+1].test(2*c+1)))
		{
			checky = 1;
		}
		else if(!boolarr[a][b+1].test(2*c+1) && (!boolarr[a][b-1].test(2*c) && boolarr[a][b-1].test(2*c+1)))
		{
			checky = -1;	
		}

		// check if (a,b,c-1) is internal(*0) and (a,b,c-1) is external(01) or viceversa
		if(!boolarr[a][b].test(2*(c-1)+1) && (!boolarr[a][b].test(2*(c+1)) && boolarr[a][b].test(2*(c+1)+1)))
		{
			checkz = 1;
		}
		else if(!boolarr[a][b].test(2*(c+1)+1) && (!boolarr[a][b].test(2*(c-1)) && boolarr[a][b].test(2*(c-1)+1)))
		{
			checkz = -1;	
		}

		//cout<<checkx<<" "<<checky<<" "<<checkz<<"\n";

		// If I didn't find the grdient along x or y direction
		if(checkx==0 || checky==0)
		{
			if(checkx==0)
				tx = 1;
			else
				tx = 0;

			if(checky==0)
				ty = 1;
			else 
				ty = 0;


			// check if (a-1,b-1,c) is internal(*0) and (a+1,b+1,c) is external(01) or viceversa
			if(!boolarr[a-1][b-1].test(2*c+1) && (!boolarr[a+1][b+1].test(2*c) && boolarr[a+1][b+1].test(2*c+1)))
			{
				checkx += tx;
				checky += ty;
			}
			else if(!boolarr[a+1][b+1].test(2*c+1) && (!boolarr[a-1][b-1].test(2*c) && boolarr[a-1][b-1].test(2*c+1)))
			{
				checkx += -tx;
				checky += -ty;
			}


			// check if (a-1,b+1,c) is internal(*0) and (a+1,b-1,c) is external(01) or viceversa
			if(!boolarr[a-1][b+1].test(2*c+1) && (!boolarr[a+1][b-1].test(2*c) && boolarr[a+1][b-1].test(2*c+1)))
			{
				checkx += tx;
				checky += -ty;
			}
			else if(!boolarr[a+1][b-1].test(2*c+1) && (!boolarr[a-1][b+1].test(2*c) && boolarr[a-1][b+1].test(2*c+1)))
			{
				checkx += -tx;
				checky += ty;
			}

		}

		// If I didn't find the grdient along x or z direction
		if(checkx==0 || checkz==0)
		{

			if(checkx==0)
				tx = 1;
			else
				tx = 0;

			if(checkz==0)
				tz = 1;
			else 
				tz = 0;



			// check if (a-1,b,c-1) is internal(*0) and (a+1,b,c+1) is external(01) or viceversa
			if(!boolarr[a-1][b].test(2*c-1) && (!boolarr[a+1][b].test(2*c+2) && boolarr[a+1][b].test(2*c+3)))
			{
				checkx += tx;
				checkz += tz;
			}
			else if(!boolarr[a+1][b].test(2*c+3) && (!boolarr[a-1][b].test(2*c-2) && boolarr[a-1][b].test(2*c-1)))
			{
				checkx += -tx;
				checkz += -tz;
			}

			// check if (a-1,b,c+1) is internal(*0) and (a+1,b,c-1) is external(01) or viceversa
			if(!boolarr[a-1][b].test(2*c+3) && (!boolarr[a+1][b].test(2*c-2) && boolarr[a+1][b].test(2*c-1)))
			{
				checkx += tx;
				checkz += -tz;
			}
			else if(!boolarr[a+1][b].test(2*c-1) && (!boolarr[a-1][b].test(2*c+2) && boolarr[a-1][b].test(2*c+3)))
			{
				checkx += -tx;
				checkz += tz;
			}
		}


		// If I didn't find the grdient along y or z direction
		if(checky==0 || checkz==0)
		{

			if(checky==0)
				ty = 1;
			else
				ty = 0;

			if(checkz==0)
				tz = 1;
			else 
				tz = 0;



			// check if (a,b-1,c-1) is internal(*0) and (a,b+1,c+1) is external(01) or viceversa
			if(!boolarr[a][b-1].test(2*c-1) && (!boolarr[a][b+1].test(2*c+2) && boolarr[a][b+1].test(2*c+3)))
			{
				checky += ty;
				checkz += tz;
			}
			else if(!boolarr[a][b+1].test(2*c+3) && (!boolarr[a][b-1].test(2*c-2) && boolarr[a][b-1].test(2*c-1)))
			{
				checky += -ty;
				checkz += -tz;
			}

			// check if (a,b-1,c+1) is internal(*0) and (a,b+1,c-1) is external(01) or viceversa
			if(!boolarr[a][b-1].test(2*c+3) && (!boolarr[a][b+1].test(2*c-2) && boolarr[a][b+1].test(2*c-1)))
			{
				checky += ty;
				checkz += -tz;
			}
			else if(!boolarr[a][b+1].test(2*c-1) && (!boolarr[a][b-1].test(2*c+2) && boolarr[a][b-1].test(2*c+3)))
			{
				checky += -ty;
				checkz += tz;
			}
		}


		// Again If I didn't find the grdient along x or y or z direction
		if(checkx==0 || checky==0 || checkz==0)
		{
			if(checkx==0)
				tx = 1;
			else
				tx = 0;

			if(checkz==0)
				tz = 1;
			else 
				tz = 0;

			if(checky==0)
				ty = 1;
			else 
				ty = 0;


			// check if (a-1,b-1,c-1) is internal(*0) and (a+1,b+1,c+1) is external(01) or viceversa
			if(!boolarr[a-1][b-1].test(2*c-1) && (!boolarr[a+1][b+1].test(2*c+2) && boolarr[a+1][b+1].test(2*c+3)))
			{
				checkx += tx;
				checky += ty;
				checkz += tz;
			}
			else if(!boolarr[a+1][b+1].test(2*c+3) && (!boolarr[a-1][b-1].test(2*c-2) && boolarr[a-1][b-1].test(2*c-1)))
			{
				checkx += -tx;
				checkz += -tz;
				checky += -ty;
			}

			// check if (a-1,b+1,c+1) is internal(*0) and (a+1,b-1,c-1) is external(01) or viceversa
			if(!boolarr[a-1][b+1].test(2*c+3) && (!boolarr[a+1][b-1].test(2*c-2) && boolarr[a+1][b-1].test(2*c-1)))
			{
				checkx += tx;
				checkz += -tz;
				checky += -ty;
			}
			else if(!boolarr[a+1][b-1].test(2*c-1) && (!boolarr[a-1][b+1].test(2*c+2) && boolarr[a-1][b+1].test(2*c+3)))
			{
				checkx += -tx;
				checkz += tz;
				checky += ty;
			}

			// check if (a+1,b-1,c+1) is internal(*0) and (a-1,b+1,c-1) is external(01) or viceversa
			if(!boolarr[a+1][b-1].test(2*c+3) && (!boolarr[a-1][b+1].test(2*c-2) && boolarr[a-1][b+1].test(2*c-1)))
			{
				checkx += -tx;
				checkz += -tz;
				checky += ty;
			}
			else if(!boolarr[a-1][b+1].test(2*c-1) && (!boolarr[a+1][b-1].test(2*c+2) && boolarr[a+1][b-1].test(2*c+3)))
			{
				checkx += tx;
				checkz += tz;
				checky += -ty;
			}

			// check if (a+1,b+1,c-1) is internal(*0) and (a-1,b-1,c+1) is external(01) or viceversa
			if(!boolarr[a+1][b+1].test(2*c-1) && (!boolarr[a-1][b-1].test(2*c+2) && boolarr[a-1][b-1].test(2*c+3)))
			{
				checkx += -tx;
				checkz += tz;
				checky += -ty;
			}
			else if(!boolarr[a-1][b-1].test(2*c+3) && (!boolarr[a+1][b+1].test(2*c-2) && boolarr[a+1][b+1].test(2*c-1)))
			{
				checkx += tx;
				checkz += -tz;
				checky += ty;
			}


		}

		//cout<<checkx<<" "<<checky<<" "<<checkz<<"\n";


		if(checkx < 0)
			check[pointarr[i].x][pointarr[i].y][pointarr[i].z].n_x *= -1;
		if(checky < 0)
			check[pointarr[i].x][pointarr[i].y][pointarr[i].z].n_y *= -1;
		if(checkz < 0)
			check[pointarr[i].x][pointarr[i].y][pointarr[i].z].n_z *= -1;

		//cout<<pointarr[i].n_x<<" "<<pointarr[i].n_y<<" "<<pointarr[i].n_z<<"\n";
		//cout<<pointarr[i].x<<" "<<pointarr[i].y<<" "<<pointarr[i].z<<"\n";
		
	}




    
	
	//x_i = 0;
	//y_i = 0;
	//z_i = 0;

	

	minframe_y = maxframe_y = y_c + (maxy - y_c)*((double)( z_i - z_c )/( maxz - z_c));
	minframe_x = maxframe_x = x_c + (maxx - x_c)*((double)( z_i - z_c )/( maxz - z_c));

	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
			for(k=0;k<2;k++)
			{
				x_i = corner + (maxx - corner)*i;
				y_i = corner + (maxy - corner)*j;
				temp_z = corner + (maxz - corner)*k;

				temp_x = x_c + (x_i - x_c)*((double)( z_i - z_c )/( temp_z - z_c));
				temp_y = y_c + (y_i - y_c)*((double)( z_i - z_c )/( temp_z - z_c)); 

				cout<<temp_x<<" "<<temp_y<<"\n";
				if(temp_x>maxframe_x)
					maxframe_x = temp_x;
				if(temp_y>maxframe_y)
					maxframe_y = temp_y;

				if(temp_x<minframe_x)
					minframe_x = temp_x;
				if(temp_y<minframe_y)
					minframe_y = temp_y;
			}

	x_left = (int)minframe_x - image_pad;
	y_left = (int)minframe_y - image_pad;

	x_right = (int)maxframe_x + image_pad;
	y_right = (int)maxframe_y + image_pad;



	cout<<"Ray casting\n";
	//bool present[400][400];

	//for(i=0;i<400;i++)
	//	for(j=0;j<400;j++)
	//		present[i][j]=false;

	// dynamically allocate memory for normal calculation
	// Rr = Ri - 2 N (Ri . N)
	present = (double **)malloc((x_right-x_left+1)*sizeof(double *));
	for(i=0;i<x_right-x_left+1;i++)
		present[i] = (double *)malloc((y_right-y_left+1)*sizeof(double));


	for(i=x_left;i<=x_right;i++)
		for(j=y_left;j<=y_right;j++)
			present[i-x_left][j-y_left] = 0.0;

	max_intens = 0.0;
	min_intens = 1000.0;

	cout<<"Input all the required variables for intensity calculation\n";
	cout<<"K_a:";
	cin>>K_a;
	cout<<"K_d:";
	cin>>K_d;
	cout<<"K_s:";
	cin>>K_s;
	cout<<"I_a:";
	cin>>I_a;
	cout<<"I_p:";
	cin>>I_p;



	// loop through all pixels in image plane
	for(i=x_left; i<=x_right; i++)
		for(j=y_left; j<=y_right; j++){

			cout<<"image pixel: "<<i<<" "<<j<<"\n";

			for(k=maxz;k>=0;k--)
			{
				
				temp_x = x_c + (double)((k - z_c)*(i - x_c))/(z_i - z_c);
				temp_y = y_c + (double)((k - z_c)*(j - y_c))/(z_i - z_c);

				//cout<< temp_x <<" "<<temp_y<<" "<<k<<"\n";
				a = (int) temp_x;
				b = (int) temp_y;
				//cout<<"point "<<a<<" "<<b<<"\n";

				bool cont = false;
				double intens = K_a * I_a, norm_x = 0.0, norm_y = 0.0, norm_z = 0.0, tot_wt = 0.0, distant;
				
				for(int no=0;no<2;no++){
					for(int kno=0;kno<2;kno++){
						if(a+no < maxsize && a+no>=0 && b+kno < maxsize && b+kno>=0  && boolarr[a+no][b+kno].test(2*k))
						{
							// check[a+no][b+kno][k] = true;
							if(dist(a+no,temp_x,b+kno,temp_y)<=1e-10)
								distant = 1e-10;
							else
								distant = dist(a+no,temp_x,b+kno,temp_y);

							norm_x += (1.0/distant)*check[a+no][b+kno][k].n_x;
							norm_y += (1.0/distant)*check[a+no][b+kno][k].n_y;
							norm_z += (1.0/distant)*check[a+no][b+kno][k].n_z;
							tot_wt += (1.0/distant);

							cont = true;
						}
					}
				}

				if(cont){
					// present[i-x_left][j-y_left] = true;
					// weighted normal
					norm_x /= tot_wt;
					norm_y /= tot_wt;
					norm_z /= tot_wt;

					normals normal(norm_x,norm_y,norm_z),incident(x_c - temp_x, y_c - temp_y, z_c - temp_z), reflected, view(x_v - temp_x, y_v - temp_y, z_v - temp_z);
					double cos_theta, cos_alpha, d_p;
					// Rr = Ri - 2 N (Ri . N)
					reflected.n_x = -incident.n_x + 2*dot_prod(incident,normal)*normal.n_x;
					reflected.n_y = -incident.n_y + 2*dot_prod(incident,normal)*normal.n_y;
					reflected.n_z = -incident.n_z + 2*dot_prod(incident,normal)*normal.n_z;

					d_p = sqrt((x_c-temp_x)*(x_c-temp_x) + (y_c-temp_y)*(y_c-temp_y) + (z_c - temp_z)*(z_c - temp_z));
					cos_theta = dot_prod(normal,incident);
					cos_alpha = dot_prod(reflected,view);

					intens += (1.0/(d_p*d_p))*I_p*(K_d*cos_theta + K_s*pow(cos_alpha,smoothness));
					present[i-x_left][j-y_left] = intens;
					if(intens > max_intens)
						max_intens = intens;
					if(intens < min_intens)
						min_intens = intens;
					break;
				}

			}
		}
	
	cout<<"camera position"<<" "<<x_c<<" "<<y_c<<" "<<z_c<<"\n";
	cout<<"Image plane: "<<x_left<<" "<<y_left<<" "<<x_right<<" "<<y_right<<" "<<z_i<<"\n";
	
	/*

	FILE *fp = fopen("first.txt","w");
	fprintf(fp, "P3\n%d %d\n255\n", x_right - x_left + 1, y_right - y_left + 1);
	for(i=x_left;i<=x_right;i++){
		for(j=y_left;j<=y_right;j++)
		{
			unsigned int color[3];
			color[0] = 0;
			color[1] = 0;
			color[2] = 0;

			if(present[i-x_left][j-y_left])
			{
				color[0] = 255;
			}
			fprintf(fp, "%u %u %u      ",color[0],color[1],color[2]);

		}
		fprintf(fp, "\n");
	}
	(void) fclose(fp);
	
	*/

	final_step(present, min_intens,max_intens,x_left,x_right,y_left,y_right);

	for(i=x_left;i<=x_right;i++)
	{
		for(j=y_left;j<=y_right;j++)
		{
			if(present[i-x_left][j-y_left]!=0.0)
				cout<<"#";
			else cout<<".";
		}
		cout<<"\n";
	}

	for(i=0;i<y_right-y_left+1;i++)
		free(present[i]);
	free(present);
	return 0;

}
