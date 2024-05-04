#include <thread>
#include <iostream>
#include <barrier>
#include <vector>
#include <bits/stdc++.h> 

using namespace std;

vector<int> getTask(int max){
    vector<int> returnVector;
    for(int index=0;index < max;index++) {
        returnVector.emplace_back(index);
    }
    return returnVector;
}

void initListOfVector(list<vector<int>> &listOfVectorTask, int numVector) {
    for(int index=0;index<numVector;index++) {
        vector<int> vector;
        listOfVectorTask.push_back(vector);
    }

    list<vector<int>>::iterator it = listOfVectorTask.begin(); 
  
    // Move the iterator by 5 elements 
    advance(it, 0); 

    vector<int> pippo;
    pippo.emplace_back(3);
    
    replace(listOfVectorTask.begin(), listOfVectorTask.end(), *it, pippo);
    advance(it, 1);
    pippo.emplace_back(4); 
    replace(listOfVectorTask.begin(), listOfVectorTask.end(), *it, pippo);
    
}

void printList(list<vector<int> >& listOfVectors) 
{ 
    for (auto vect : listOfVectors) { 
        // Each element of the list is 
        // a vector itself 
        vector<int> currentVector = vect; 
  
        cout << "[ "; 
  
        // Printing vector contents 
        for (auto element : currentVector) 
            cout << element << ' '; 
  
        cout << ']'; 
        cout << '\n'; 
    } 
} 

int main(int argc, char const *argv[])
{
    list<vector<int>> listOfVectorTask;
    vector<thread> threads;
    int numThreads = 3;
    
    initListOfVector(listOfVectorTask, numThreads);
    
    printList(listOfVectorTask);
    
    for(int index=0;index<numThreads;index++) {

    }
}
