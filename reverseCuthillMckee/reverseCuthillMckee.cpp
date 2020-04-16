// reverseCuthillMckee.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include<fstream>
#include <vector>
#include <string>

using namespace std;

void matrixAdjasence(int &num_adj, std::vector<int> &xadj, std::vector<int> &adj, std::vector<std::vector<bool>> A)
{
	const int num_node = A.size();
	std::vector<std::vector<bool>> AwithoutDiag(num_node, std::vector<bool>(num_node, 0));
	for (int i = 0; i < A.size(); ++i)
	{
		for (int j = 0; j < num_node; ++j)
		{
			if (i != j)
			{
				AwithoutDiag[i][j] = A[i][j];
			}
		}
	}
	for (int i = 0; i < num_node; ++i)
	{
		std::vector<int> indexes;
		for (int j = 0; j < num_node; ++j)
		{
			if (AwithoutDiag[i][j] != 0)
			{
				indexes.push_back(j);
			}
		}
		xadj.push_back(adj.size());
		for (int k = 0; k < indexes.size(); ++k)
		{
			adj.push_back(indexes[k]);
		}
	}
	xadj.push_back(adj.size());
	num_adj = adj.size();
}

void rootls(int root, int num_adj, std::vector<int> xadj, std::vector<int> adj,
	int num_node, std::vector<bool> &mask, int &num_ls, std::vector<int> &xls, std::vector<int> &ls)
{
	// Define ls structure with root in ROOT.
	mask[root] = false;
	xls.clear();
	ls.clear();
	ls.push_back(root);
	num_ls = -1;
	int l_end = -1;
	int num_node_ls = 0;

	// l_begin is a pointer to the beginning of the current level, 
	// and l_end indicates the end of this level

	while (true)
	{
		int l_begin = l_end + 1;
		l_end = num_node_ls;
		num_ls++;
		xls.push_back(l_begin);

		// Generate the next level by finding all the masked neighbors of nodes
		// in the current level

		for (int i = l_begin; i <= l_end; ++i)
		{
			int node = ls[i];
			int jstrt = xadj[node];
			int jstop = xadj[node + 1] - 1;
			for (int j = jstrt; j <= jstop; ++j)
			{
				int voisin = adj[j];

				if (mask[voisin])
				{
					num_node_ls++;
					ls.push_back(voisin);
					mask[voisin] = 0;
				}
			}
		}


		int l_size = num_node_ls - l_end;
		if (l_size <= 0)
		{
			break;
		}
	}

	xls.push_back(l_end + 1);
	for (int i = 0; i < ls.size(); ++i)
	{
		mask[ls[i]] = true;
	}
}

void fnroot(int num_adj, std::vector<int> xadj, std::vector<int> adj,
	std::vector<bool> mask, int num_node, int &root, int &num_ls, std::vector<int> &xls, std::vector<int> &ls)
{
	int num_ls_2;
	rootls(root, num_adj, xadj, adj, num_node,
		mask, num_ls, xls, ls);

	int num_node_ls = xls[num_ls + 1] - 1;

	if (num_ls == 0)
	{
		return;
	}
	if (num_ls == num_node_ls)
	{
		return;
	}
	while (true)
	{
		int mindeg = num_node_ls + 1;
		int jstrt = xls[num_ls];
		root = ls[jstrt];

		if (jstrt < num_node_ls)
		{
			for (int j = jstrt; j <= num_node_ls; ++j)
			{
				int node = ls[j];
				int ndeg = 0;
				int kstrt = xadj[node];
				int kstop = xadj[node + 1] - 1;

				for (int k = kstrt; k <= kstop; ++k)
				{
					int voisin = adj[k];
					if (mask[voisin])
					{
						ndeg++;
					}
				}
				if (ndeg < mindeg)
				{
					root = node;
					mindeg = ndeg;
				}
			}
		}


		rootls(root, num_adj, xadj, adj, num_node,
			mask, num_ls, xls, ls);

		num_node_ls = xls[num_ls + 1] - 1;

		if (num_ls == 0)
		{
			return;
		}
		if (num_ls == num_node_ls)
		{
			return;
		}

		while (true)
		{
			int mindeg = num_node_ls + 1;
			int jstrt = xls[num_ls];
			root = ls[jstrt];

			if (jstrt < num_node_ls)
			{
				for (int j = jstrt; j <= num_node_ls; ++j)
				{
					int node = ls[j];
					int ndeg = 0;
					int kstrt = xadj[node];
					int kstop = xadj[node + 1] - 1;

					for (int k = kstrt; k <= kstop; ++k)
					{
						int voisin = adj[k];
						if (mask[voisin])
						{
							ndeg++;
						}
					}
					if (ndeg < mindeg)
					{
						root = node;
						mindeg = ndeg;
					}
				}
			}

			
			rootls(root, num_adj, xadj, adj, num_node,
				mask, num_ls_2, xls, ls);

			if (num_ls_2 <= num_ls)
			{
				break;
			}
			num_ls = num_ls_2;

			if (num_node_ls <= num_ls)
			{
				break;
			}
		}

		rootls(root, num_adj, xadj, adj, num_node,
			mask, num_ls_2, xls, ls);
		if (num_ls_2 <= num_ls)
		{
			break;
		}
		num_ls = num_ls_2;
		if (num_node_ls <= num_ls)
		{
			break;
		}

	}

	

}

void degree(int root, int num_adj, std::vector<int> xadj, std::vector<int> adj, std::vector<bool> mask, int num_node,
	std::vector<int> &deg, int &num_node_ls, std::vector<int> &ls)
{
	// The sign of xadj(I) is used to indicate if node I has been considered.
	ls.clear();
	ls.push_back(root);
	xadj[root] = -xadj[root];
	int lvlend = -1;
	int node;
	num_node_ls = 0;

	// LBEGIN is the pointer to the beginning of the current ls, and
	// LVLEND points to the end of this ls.

	while (true)
	{
		int lbegin = lvlend + 1;
		lvlend = num_node_ls;

		// Find the degrees of nodes in the current ls,
		// and at the same time, generate the next ls.

		for (int i = lbegin; i <= lvlend; ++i)
		{
			node = ls[i];
			int jstrt = -xadj[node];
			int jstop = abs(xadj[node + 1]) - 1;
			int ideg = 0;

			for (int j = jstrt; j <= jstop; ++j)
			{
				int nbr = adj[j];
				if (mask[nbr])
				{
					ideg++;
					if (0 < xadj[nbr])
					{
						xadj[nbr] = -xadj[nbr];
						num_node_ls++;
						ls.push_back(nbr);
					}
				}
			}
			deg.resize(node + 1);
			deg[node] = ideg;
		}

		// Compute the current ls width

		int lvsize = num_node_ls - lvlend;

		// If the current ls width is nonzero, generate another ls

		if (lvsize == 0)
		{
			break;
		}	
	}

	// Reset xadj to its correct sign and return

	for (int i = 0; i <= num_node_ls; ++i)
	{
		node = ls[i];
		xadj[node] = -xadj[node];
	}

}

void rcm(int root, int num_adj, std::vector<int> xadj, std::vector<int> adj, int num_node,
	std::vector<bool> &mask, std::vector<int> &perm, int &num_node_ls)
{
	// Make sure num_node is legal.
	if (num_node < 1)
	{
		cout << endl << "RCM - Fatal error!" << endl <<
			"Illegal input value of num_node = " << num_node << endl <<
			"Acceptable values must be positive." << endl;
		exit(1);
	}
	// Make sure ROOT is legal.

	if (root < 0 || num_node < root)
	{
		cout << endl << "RCM - Fatal error!" << endl <<
			"Illegal input value of ROOT = " << root << endl <<
			"Acceptable values are between 0 and " << num_node - 1 << endl;
		exit(1);
	}

	// Find the degrees of the nodes in the component specified by MASK and ROOT.

	std::vector<int> deg;
	degree(root, num_adj, xadj, adj, mask, num_node,
		deg, num_node_ls, perm);

	mask[root] = false;

	if (num_node_ls < 0)
	{
		cout << endl << "RCM - Fatal error!" << endl <<
			"Inexplicable component size num_node_ls = " << num_node_ls << endl;
		exit(1);
	}

	// If the connected component is a singleton, there is no reordering to do.
	if (num_node_ls == 0)
	{
		return;
	}

	int lvlend = -1;
	int lnbr = 0;

	while (lvlend < lnbr)
	{
		int lbegin = lvlend + 1;
		lvlend = lnbr;

		for (int i = lbegin; i <= lvlend; ++i)
		{
			int node = perm[i];
			int jstrt = xadj[node];
			int jstop = xadj[node + 1] - 1;

			// FNBR and LNBR point to the first and last neighbors
			// of the current node in PERM.

			int fnbr = lnbr + 1;

			for (int j = jstrt; j <= jstop; ++j)
			{

				int nbr = adj[j];
				if (mask[nbr])
				{
					lnbr++;
					mask[nbr] = false;
					perm[lnbr] = nbr;
				}
			}

			// If no neighbors, skip to next node in this ls.
			if (lnbr <= fnbr)
			{
				continue;
			}

			// Sort the neighbors of NODE in increasing order by degree.
			// Linear insertion is used.

			int k = fnbr;

			while (k < lnbr)
			{
				int l = k;
				k++;
				int nbr = perm[k];

				while (fnbr < l)
				{
					int lperm = perm[l];

					if (deg[lperm] <= deg[nbr])
					{
						break;
					}
					perm[l + 1] = lperm;
					l = l - 1;
				}
				perm[l + 1] = nbr;
			}
		}
		
	}

	// We now have the Cuthill - McKee ordering.
	// Reverse it to get the Reverse Cuthill - McKee ordering.
	
	std::vector <int> temp;
	for (int i = 0; i < perm.size(); ++i)
	{
		temp.push_back(perm[perm.size() - i - 1]);
	}
	perm = temp;
		
}

void symrcmhw(std::vector<std::vector<bool>> A, std::vector<int>&perm)
{
	perm.clear();
	const int num_node = A.size();
	int num_adj;
	int num = 0;
	int root;
	std::vector<int> xadj;
	std::vector<int> adj;
	std::vector<bool> mask(num_node, true);

	matrixAdjasence(num_adj, xadj, adj, A);

	for (int i = 0; i < num_node; ++i)
	{
		if (mask[i])
		{
			root = i;
			// Find a pseudo - peripheral node ROOT.The level structure found by
			// ROOT_FIND is stored starting at PERM(NUM).
			int num_ls;
			std::vector<int> xls;
			std::vector<int> ls;
			fnroot(num_adj, xadj, adj,
				mask, num_node, root, num_ls, xls, ls);

			// RCM orders the component using ROOT as the starting node.
			int num_node_ls;
			rcm(root, num_adj, xadj, adj, num_node,
				mask, ls, num_node_ls);

			perm.resize(perm.size()+ls.size());
			for (int i = num; i < num + num_node_ls + 1; ++i)
			{
				perm[i] = ls[i - num];
			}

			num = num + num_node_ls + 1;
			if (num_node <= num)
			{
				return;
			}
		}
	}
}

void readTriplets(std::vector<std::vector<float>> &triplets, string fileName)
{
	int n, m, num;

	std::ifstream infile;
	infile.open(fileName);
	infile >> n;
	infile >> m;
	infile >> num;

	triplets.resize(num + 1);
	for (int i = 0; i < triplets.size(); ++i)
	{
		triplets[i].resize(3);
	}

	triplets[0][0] = n;
	triplets[0][1] = m;
	triplets[0][2] = num;

	int idx = 1;
	while (!infile.eof()) 
	{
		int x, y, val;
		infile >> x;
		infile >> y;
		infile >> val;

		triplets[idx][0] = x;
		triplets[idx][1] = y;
		triplets[idx][2] = val;
		idx++;
	}
	infile.close();
}

void triplets2matrix(std::vector<std::vector<float>> triplets, std::vector<std::vector<bool>> &mat)
{
	int n = triplets[0][0];
	int m = triplets[0][1];
	int num = triplets[0][2];

	if (n != m)
	{
		cout << "The matrix must be square" << endl;
		exit(1);
	}

	mat.resize(n);
	for (int i = 0; i < n; ++i)
	{
		mat[i].resize(m);
	}

	for (int i = 1; i < num + 1; ++i)
	{
		int x = triplets[i][0];
		int y = triplets[i][1];
		int val = triplets[i][2];

		mat[x-1][y-1] = val;
	}
}

void writePermutation(std::vector<int> perm, string fileName)
{
	std::ofstream outfile(fileName, ios::trunc);

	for (int i = 0; i < perm.size(); ++i)
	{
		outfile << perm[i] << endl;
	}
	outfile.close();
}

int main()
{
	const string fileNameIn = "bucky.txt";
	const string fileNameOut = "bucky_perm.txt";
	std::vector<std::vector<float>> triplets;
	std::vector<std::vector<bool>> mat;
	std::vector<int> perm;

	readTriplets(triplets, fileNameIn);
	triplets2matrix(triplets, mat);

	
	symrcmhw(mat, perm);

	writePermutation(perm, fileNameOut);
	
	return 0;
}
