#include <H5Apublic.h>
#include <H5Ipublic.h>
#include <H5Tpublic.h>
#include <stdlib.h>
#include <hdf5.h>
#include "auxi.h"
#include "struct.h"

#define TRIG_D_NAME   "triggers"
#define RANK          1

int trig_h5_init (char *outname, Command_line_opts *opts, Search_settings *sett, 
     Search_range *s_range, Trigger *sgnlv)
{
     int        i, j;
     hid_t      file, t_dataset, t_space, t_prop, filespace, memspace, space_scalar;
     hsize_t    t_dim[RANK] = {0};
     hsize_t    t_maxdim[RANK] = {H5S_UNLIMITED};
     herr_t     status;
     
     hid_t vstr_type_id = H5Tcopy(H5T_C_S1);
     H5Tset_size(vstr_type_id, H5T_VARIABLE);
     
     // ------------------------------------------------------------------------
     // Define data types
     
     // triggers data type (sgnlv)
     hid_t ffstat_type = H5Tvlen_create(H5T_NATIVE_FLOAT);
     hid_t t_tid = H5Tcreate (H5T_COMPOUND, sizeof(Trigger));
     H5Tinsert(t_tid, "m", HOFFSET(Trigger, m), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "n", HOFFSET(Trigger, n), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "s", HOFFSET(Trigger, s), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "ra", HOFFSET(Trigger, ra), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "dec", HOFFSET(Trigger, dec), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "fdot", HOFFSET(Trigger, fdot), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "ffstat", HOFFSET(Trigger, ffstat), ffstat_type);

     // Command line options data type  (opts)
     hid_t cmd_opts_tid = H5Tcreate (H5T_COMPOUND, sizeof(Command_line_opts));
     H5Tinsert(cmd_opts_tid, "checkp_flag", HOFFSET(Command_line_opts, checkp_flag), H5T_NATIVE_INT);
     H5Tinsert(cmd_opts_tid, "veto_flag", HOFFSET(Command_line_opts, veto_flag), H5T_NATIVE_INT);
     //H5Tinsert(cmd_opts_tid, "gen_vlines_flag", HOFFSET(Command_line_opts, gen_vlines_flag), H5T_NATIVE_INT);
     //H5Tinsert(cmd_opts_tid, "help_flag", HOFFSET(Command_line_opts, help_flag), H5T_NATIVE_INT);
     H5Tinsert(cmd_opts_tid, "seg", HOFFSET(Command_line_opts, seg), H5T_NATIVE_INT);
     H5Tinsert(cmd_opts_tid, "band", HOFFSET(Command_line_opts, band), H5T_NATIVE_INT);
     H5Tinsert(cmd_opts_tid, "hemi", HOFFSET(Command_line_opts, hemi), H5T_NATIVE_INT);
     H5Tinsert(cmd_opts_tid, "nod", HOFFSET(Command_line_opts, nod), H5T_NATIVE_INT);
     H5Tinsert(cmd_opts_tid, "thr", HOFFSET(Command_line_opts, thr), H5T_NATIVE_DOUBLE);
     //H5Tinsert(cmd_opts_tid, "fpo_val", HOFFSET(Command_line_opts, fpo_val), H5T_NATIVE_DOUBLE);
     H5Tinsert(cmd_opts_tid, "narrowdown", HOFFSET(Command_line_opts, narrowdown), H5T_NATIVE_DOUBLE);
     H5Tinsert(cmd_opts_tid, "overlap", HOFFSET(Command_line_opts, overlap), H5T_NATIVE_DOUBLE);
     H5Tinsert(cmd_opts_tid, "indir", HOFFSET(Command_line_opts, indir), vstr_type_id);
     H5Tinsert(cmd_opts_tid, "outdir", HOFFSET(Command_line_opts, outdir), vstr_type_id);
     H5Tinsert(cmd_opts_tid, "range_file", HOFFSET(Command_line_opts, range_file), vstr_type_id);
     H5Tinsert(cmd_opts_tid, "grid_file", HOFFSET(Command_line_opts, grid_file), vstr_type_id);
     //H5Tinsert(cmd_opts_tid, "dump_range_file", HOFFSET(Command_line_opts, dump_range_file), H5T_C_S1);
     H5Tinsert(cmd_opts_tid, "usedet", HOFFSET(Command_line_opts, usedet), vstr_type_id);
     H5Tinsert(cmd_opts_tid, "addsig", HOFFSET(Command_line_opts, addsig), vstr_type_id);
     H5Tinsert(cmd_opts_tid, "fstat_norm", HOFFSET(Command_line_opts, fstat_norm), vstr_type_id);
     //H5Tinsert(cmd_opts_tid, "label", HOFFSET(Command_line_opts, label), H5T_C_S1);
     //H5Tinsert(cmd_opts_tid, "state_file", HOFFSET(Command_line_opts, state_file), H5T_C_S1);

     // Search settings data type  (sett)
     
     printf("sett size=%d\n", sizeof(Search_settings)-sizeof(sett->lines)+
          sett->numlines_band*2*sizeof(double));
     hid_t sett_tid = H5Tcreate (H5T_COMPOUND, sizeof(Search_settings)-sizeof(sett->lines)+
          sett->numlines_band*2*sizeof(double));
     H5Tinsert(sett_tid, "fpo", HOFFSET(Search_settings, fpo), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "dt", HOFFSET(Search_settings, dt), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "B", HOFFSET(Search_settings, B), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "oms", HOFFSET(Search_settings, oms), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "omr", HOFFSET(Search_settings, omr), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "Smin", HOFFSET(Search_settings, Smin), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "Smax", HOFFSET(Search_settings, Smax), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "sepsm", HOFFSET(Search_settings, sepsm), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "cepsm", HOFFSET(Search_settings, cepsm), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "nfft", HOFFSET(Search_settings, nfft), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nod", HOFFSET(Search_settings, nod), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "N", HOFFSET(Search_settings, N), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nfftf", HOFFSET(Search_settings, nfftf), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nmax", HOFFSET(Search_settings, nmax), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nmin", HOFFSET(Search_settings, nmin), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "s", HOFFSET(Search_settings, s), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nd", HOFFSET(Search_settings, nd), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "interpftpad", HOFFSET(Search_settings, interpftpad), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "fftpad", HOFFSET(Search_settings, fftpad), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "Ninterp", HOFFSET(Search_settings, Ninterp), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nifo", HOFFSET(Search_settings, nifo), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "numlines_band", HOFFSET(Search_settings, numlines_band), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "bufsize", HOFFSET(Search_settings, bufsize), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "dd", HOFFSET(Search_settings, dd), H5T_NATIVE_INT);
     hsize_t dims[2] = {sett->numlines_band, 2};
     hid_t lines_t = H5Tarray_create2(H5T_NATIVE_DOUBLE, 2, dims);
     H5Tinsert(sett_tid, "lines", HOFFSET(Search_settings, lines), lines_t);
     
     // ------------------------------------------------------------------------

     
     //Create the file.
     file = H5Fcreate(outname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
     
     hid_t scalar_space_id = H5Screate(H5S_SCALAR);
     hid_t attr = H5Acreate2(file, "opts", cmd_opts_tid, scalar_space_id, H5P_DEFAULT, H5P_DEFAULT);
     status = H5Awrite(attr, cmd_opts_tid, opts);
     H5Aclose(attr);
     
     attr = H5Acreate2(file, "sett", sett_tid, scalar_space_id, H5P_DEFAULT, H5P_DEFAULT);
     status = H5Awrite(attr, sett_tid, sett);
     H5Aclose(attr);

     hsize_t chunk_dims[1] = {3};
     t_prop = H5Pcreate(H5P_DATASET_CREATE);
     status = H5Pset_chunk(t_prop, RANK, chunk_dims);
     
     // Create the data space.
     t_space = H5Screate_simple(RANK, t_dim, t_maxdim);
     
     // Create the dataset.
     t_dataset = H5Dcreate2(file, TRIG_D_NAME, t_tid, t_space, H5P_DEFAULT, t_prop, H5P_DEFAULT);
     
     // Wtite data to the dataset;
     status = H5Dwrite(t_dataset, t_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, sgnlv);
     
     // Release resources
     //H5Tclose(tr_tid);
     H5Pclose(t_prop);
     H5Dclose(t_dataset);
     H5Sclose(t_space);
     H5Fclose(file);
     H5Tclose(ffstat_type);

     return(EXIT_SUCCESS);
}



int trig_h5_extend(char *outname, Trigger *sgnlv, int sgnlv_size){
     
     // Extend the triggers dataset by sgnlv_size entries from sgnlv bffer
     // Assume the file is already created and contains dataset "triggers"
     // The function closes all HDF5 objects it opens

     int        i, j;
     hid_t      file, t_dataset, t_space, t_prop, filespace, memspace, space_scalar;
     hsize_t    t_dim[RANK] = {0};
     hsize_t    t_maxdim[RANK] = {H5S_UNLIMITED};
     herr_t     status;
     
     // ------------------------------------------------------------------------
     // Define data types
     
     // triggers data type (sgnlv)
     hid_t ffstat_type = H5Tvlen_create(H5T_NATIVE_FLOAT);
     hid_t t_tid = H5Tcreate (H5T_COMPOUND, sizeof(Trigger));
     H5Tinsert(t_tid, "m", HOFFSET(Trigger, m), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "n", HOFFSET(Trigger, n), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "s", HOFFSET(Trigger, s), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "ra", HOFFSET(Trigger, ra), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "dec", HOFFSET(Trigger, dec), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "fdot", HOFFSET(Trigger, fdot), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "ffstat", HOFFSET(Trigger, ffstat), ffstat_type);

     // ------------------------------------------------------------------------

     //Open the file.
     file = H5Fopen(outname, H5F_ACC_RDWR, H5P_DEFAULT);
     
     // Open the dataset.
     t_dataset = H5Dopen2(file, TRIG_D_NAME, H5P_DEFAULT);
     
     // Get the dataspace.
     filespace = H5Dget_space(t_dataset);
     
     // Get the current size of the dataset.
     status = H5Sget_simple_extent_dims(filespace, t_dim, NULL);
     
     // Extend the dataset.
     t_dim[0] += sgnlv_size;
     status = H5Dset_extent(t_dataset, t_dim);
     
     // Close the old dataspace and get the new one.
     status = H5Sclose(filespace);
     filespace = H5Dget_space(t_dataset);
     
     // Define hyperslab in the extended portion of the dataset.
     hsize_t start[RANK] = {t_dim[0]-sgnlv_size};
     hsize_t count[RANK] = {sgnlv_size};
     status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);
     memspace = H5Screate_simple(RANK, count, NULL);
     status = H5Dwrite(t_dataset, t_tid, memspace, filespace, H5P_DEFAULT, sgnlv);
     
     // Release resources
     //H5Tclose(tr_tid);
     H5Sclose(memspace);
     H5Sclose(filespace);
     H5Tclose(ffstat_type);
     H5Dclose(t_dataset);
     H5Fclose(file);
     
     return(1);
}
