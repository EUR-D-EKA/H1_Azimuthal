#include <iostream>

int MakeQueue(){
	int nFiles, nComputers, nDigits;
	std::string name, fName, fPath;

	std::cout<<"Enter name"<<std::endl;
	std::cin >> name;

	std::cout<<"Enter number of files to process"<<std::endl;
	std::cin >> nFiles;

	std::cout<<"Enter number of computers to request"<<std::endl;
	std::cin >> nComputers;

	std::cout<<"Enter path to files"<<std::endl;
	std::cin >> fPath;

	std::cout<<"Enter file name"<<std::endl;
	std::cin >> fName;

	std::cout<<"Enter number of digits in file name"<<std::endl;
	std::cin >> nDigits;


	int step = static_cast<int>(nFiles/nComputers)-1;

	ofstream Queue("Queue.txt");
	int n = 0;
  	for (int i = 0; i < nComputers; ++i){
  		char x[1000];
  		std::string s = name + "_" + std::to_string(i);
  		sprintf(x, "../../input/%s.steer %s out_data",s.c_str(),s.c_str());
  		Queue << x <<std::endl;

  		ofstream Steerfile;
  		Steerfile.open("../input/"+s+".steer");

  		(i==0)? n = 0 : ++n;
  		std::string digit1(nDigits, '0');
  		digit1 += std::to_string( n );
  		digit1.erase(0,digit1.size() - nDigits);

  		n += step;
  		if (i == nComputers-1 && n<nFiles-1) n = nFiles-1 ;
  		std::string digit2(nDigits, '0');
  		digit2 += std::to_string( n );
  		digit2.erase(0,digit2.size() - nDigits);

  		std::string fNames = fName + digit1 + "-" + digit2;
  		sprintf(x, "H1SteerTree() {\nfPath = \"%s\";\nfFileNames= \"%s\";\nfLoadHAT = true;\nfLoadMODS = true;\nfHatSelection=\"(fRunType>0)||(fIl1ac[67]>0)\";\n}",fPath.c_str(),fNames.c_str());
  		Steerfile<<x<<std::endl;
  		Steerfile.close();
  	}
  	Queue.close();
  	return 0;
}
