void parse_dString(char* buffer, double d_out[], int n_size)
//=Description=================
//
// parse an double value line into the return array
// of size n_size
//IN   : buffer[256]
//OUT  : d_out[n_size]
//
//HEADER <strstream>
//=============================
{
  istrstream one_line(buffer); // construct an istringstream class object with buffer
                               // note: new version= istringstream
  for(int i=0; i<n_size; i++ )
  {
     d_out[i]=0; //initialize
     one_line >> d_out[i];     // one double value into the out array
  }
}


void parse_iString(char* buffer, int i_out[], int n_size)
//=Description=================
//
// parse an integer line into the return array
// of size n_size
//IN   : buffer[256]
//OUT  : i_out[n_size]
//
//HEADER <strstream>
//=============================
{
  istrstream one_line(buffer); // construct an istringstream class object with buffer
                               // note: new version= istringstream
  for(int i=0; i<n_size; i++ )
  {
     i_out[i]=0; //initialize
     one_line >> i_out[i];     // one integer value into the out array
  }
}
