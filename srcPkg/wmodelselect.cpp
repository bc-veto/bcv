#include <cstdlib>

#include <iostream>

using namespace std;

// a set of intervals with associated values (a, b) -> c

struct interval
{
    double v; // value of the interval
    int left; // if non-zero, index of the left interval contributing to the value
    int right; // if left is non-zero, index of the right interval contributing to the value
    interval() : v(0), left(0), right(0) {} // paranoid zeroing
};

int n; // the number of unique endpoints in the problem (~hundreds)
interval* p; // pointer to an [n*n] array storing intervals

// comnvert between 1D and 2D indexing into the p array
int ifromab(int a, int b)
{
    return a + n * b;
}

int afromi(int i)
{
    return i % n;
}

int bfromi(int i)
{
    return i / n;
}

// recursively work out what the contributors to an interval are, after the
// working array has been processed

void decompose(int i)
{
    if (p[i].left > 0)
    {
        // we need to chase
        decompose(p[i].left);
        decompose(p[i].right);
    }
    else
    {
        if (p[i].v > 0)
        {
            // this is a selected seeded interval
            cout << "(" << afromi(i) << ", " << bfromi(i) << ")" << " -> " << p[i].v << endl;
        }
        else
        {
            // no seeded interval covers this region
        }
    }
}

int main(int argc, char* argv[])
{    
    // the number of unique endpoints, mapped onto the integers
    
    n = 512;
    
    // a square array of working memory
    
    p = new interval[n * n];
    
    // seed with a few intervals
    
    /*
    p[ifromab(10, 20)].v = 1;
    p[ifromab(20, 30)].v = 1;
    p[ifromab(15, 25)].v = 1;
    p[ifromab(10, 30)].v = 1;
    p[ifromab( 5, 35)].v = 1;
    */
    
    for (int z = 0; z != 1000; ++z)
    {
        p[ifromab(rand() % n, rand() % n)].v = double(rand()) / RAND_MAX;
    }
    
    // find the set of non-overlapping intervals with the greatest total 
    // value
    
    // loop over longer and longer intervals
    for (int w = 2; w < n; ++w)
    {
        // loop over intervals with fixed length and different starting positions
        for (int a = 0; a + w < n; ++a)
        {
            int b = a + w;
                     
            int i = ifromab(a,b);
            
            // loop over midpoints partitioning into left and right subpintervals
            for (int c = a + 1; c != b; ++c)
            {
                int j = ifromab(a, c);
                int k = ifromab(c, b);
                // does the combination of the left and right intervals outperform the
                // seeded value, if any, or the previous best left and right intervals?
                if (p[i].v < p[j].v + p[k].v)
                {
                    // if so, set up the interval to use the two subintervals
                    p[i].v = p[j].v + p[k].v;
                    p[i].left = j;
                    p[i].right = k;
                }
            }            
        }
    }

    // now navigate recursively to print out the intervals
    cout << "(0, " << n - 1 << ") -> " << p[ifromab(0,n-1)].v << endl << "    from" << endl;
    decompose(ifromab(0,n-1));    
}