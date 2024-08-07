/*
 * @file cs241.c
 * @brief Main file for CS241 library.
 * @author Mei-Hui Su - SCEC <mei@usc.edu>
 * @author Brad -- update to PROJ
 * @version 1.0
 *
 * @section DESCRIPTION
 *
 * Delivers the CS241 model 
 * 
 *
 */

#include "ucvm_model_dtypes.h"
#include "cs241.h"
#include <assert.h>

/** The config of the model */
char *cs241_config_string=NULL;
int cs241_config_sz=0;

/**
 * Initializes the CS241 plugin model within the UCVM framework. In order to initialize
 * the model, we must provide the UCVM install path and optionally a place in memory
 * where the model already exists.
 *
 * @param dir The directory in which UCVM has been installed.
 * @param label A unique identifier for the velocity model.
 * @return Success or failure, if initialization was successful.
 */
int cs241_init(const char *dir, const char *label) {
    char cs241_projstr[64];
    int tempVal = 0;
    char configbuf[512];
    double north_height_m = 0, east_width_m = 0, rotation_angle = 0;

    // Initialize variables.
    cs241_configuration = calloc(1, sizeof(cs241_configuration_t));
    cs241_velocity_model = calloc(1, sizeof(cs241_model_t));
    cs241_vs30_map = calloc(1, sizeof(cs241_vs30_map_config_t));

    cs241_config_string = calloc(CS241_CONFIG_MAX, sizeof(char));
    cs241_config_string[0]='\0';
    cs241_config_sz=0;

    // Configuration file location.
    sprintf(configbuf, "%s/model/%s/data/config", dir, label);

    // Set up model directories.
    sprintf(cs241_vs30_etree_file, "%s/model/ucvm/ucvm.e", dir);

    // Read the cs241_configuration file.
    if (cs241_read_configuration(configbuf, cs241_configuration) != SUCCESS)
        return FAIL;

    // Set up the iteration directory.
    sprintf(cs241_iteration_directory, "%s/model/%s/data/%s/", dir, label, cs241_configuration->model_dir);

    // Can we allocate the model, or parts of it, to memory. If so, we do.
    tempVal = cs241_try_reading_model(cs241_velocity_model);

    if (tempVal == SUCCESS) {
        fprintf(stderr, "WARNING: Could not load model into memory. Reading the model from the\n");
        fprintf(stderr, "hard disk may result in slow performance.");
    } else if (tempVal == FAIL) {
        cs241_print_error("No model file was found to read from.");
        return FAIL;
    }

    if (cs241_read_vs30_map(cs241_vs30_etree_file, cs241_vs30_map) != SUCCESS) {
        cs241_print_error("Could not read the Vs30 map data from UCVM.");
        return FAIL;
    }

//reference ????
    char* pstr= "+proj=utm +zone=10 +ellps=clrk66 +datum=NAD27 +units=m +no_defs";

    // We need to convert the point from lat, lon to UTM, let's set it up.
    /* Setup projection */

    // We need to convert the point from lat, lon to UTM, let's set it up.
    snprintf(cs241_projstr, 64, "+proj=utm +zone=%d +datum=NAD27 +units=m +no_defs", cs241_configuration->utm_zone);
    if (!(cs241_geo2utm = proj_create_crs_to_crs(PJ_DEFAULT_CTX, "EPSG:4326", cs241_projstr, NULL))) {
        cs241_print_error("Could not set up Proj transformation from EPSG:4325 to UTM.");
        cs241_print_error((char  *)proj_context_errno_string(PJ_DEFAULT_CTX, proj_context_errno(PJ_DEFAULT_CTX)));
        return (UCVM_CODE_ERROR);
    }

    assert(cs241_vs30_map);
    if (!(cs241_geo2aeqd = proj_create_crs_to_crs(PJ_DEFAULT_CTX, "EPSG:4326", cs241_vs30_map->projection, NULL))) {
        cs241_print_error("Could not set up Proj transformation from EPSG:4326 to AEQD projection.");
        cs241_print_error((char  *)proj_context_errno_string(PJ_DEFAULT_CTX, proj_context_errno(PJ_DEFAULT_CTX)));
        return (UCVM_CODE_ERROR);
    }

    // In order to simplify our calculations in the query, we want to rotate the box so that the bottom-left
    // corner is at (0m,0m). Our box's height is cs241_total_height_m and cs241_total_width_m. We then rotate the
    // point so that is is somewhere between (0,0) and (cs241_total_width_m, cs241_total_height_m). How far along
    // the X and Y axis determines which grid points we use for the interpolation routine.

    // Calculate the rotation angle of the box.
    assert(cs241_configuration);
    north_height_m = cs241_configuration->top_left_corner_n - cs241_configuration->bottom_left_corner_n;
    east_width_m = cs241_configuration->top_left_corner_e - cs241_configuration->bottom_left_corner_e;

    // Rotation angle. Cos, sin, and tan are expensive computationally, so calculate once.
    rotation_angle = atan(east_width_m / north_height_m);

    cs241_cos_rotation_angle = cos(rotation_angle);
    cs241_sin_rotation_angle = sin(rotation_angle);

    cs241_total_height_m = sqrt(pow(cs241_configuration->top_left_corner_n - cs241_configuration->bottom_left_corner_n, 2.0f) +
                              pow(cs241_configuration->top_left_corner_e - cs241_configuration->bottom_left_corner_e, 2.0f));
    cs241_total_width_m = sqrt(pow(cs241_configuration->top_right_corner_n - cs241_configuration->top_left_corner_n, 2.0f) +
                             pow(cs241_configuration->top_right_corner_e - cs241_configuration->top_left_corner_e, 2.0f));

    // Get the cos and sin for the Vs30 map rotation.
    cs241_cos_vs30_rotation_angle = cos(cs241_vs30_map->rotation * DEG_TO_RAD);
    cs241_sin_vs30_rotation_angle = sin(cs241_vs30_map->rotation * DEG_TO_RAD);

    /* setup config_string */
    sprintf(cs241_config_string,"config = %s\n",configbuf);
    cs241_config_sz=1;

    // Let everyone know that we are initialized and ready for business.
    cs241_is_initialized = 1;

    return SUCCESS;
}

/**
 * Queries CS241 at the given points and returns the data that it finds.
 *
 * @param points The points at which the queries will be made.
 * @param data The data that will be returned (Vp, Vs, density, Qs, and/or Qp).
 * @param numpoints The total number of points to query.
 * @return SUCCESS or FAIL.
 */
int cs241_query(cs241_point_t *points, cs241_properties_t *data, int numpoints) {
    int i = 0;

    double point_u = 0, point_v = 0;
    double point_x = 0, point_y = 0;

    int load_x_coord = 0, load_y_coord = 0, load_z_coord = 0;
    double x_percent = 0, y_percent = 0, z_percent = 0;
    cs241_properties_t surrounding_points[8];

    int zone = 10;
    int longlat2utm = 0;

    for (i = 0; i < numpoints; i++) {

        // We need to be below the surface to service this query.
        if (points[i].depth < 0) {
            data[i].vp = -1;
            data[i].vs = -1;
            data[i].rho = -1;
            data[i].qp = -1;
            data[i].qs = -1;
            continue;
        }

        PJ_COORD xyzSrc = proj_coord(points[i].latitude, points[i].longitude, 0.0, HUGE_VAL);
        PJ_COORD xyzDest = proj_trans(cs241_geo2utm, PJ_FWD, xyzSrc);
        int err = proj_context_errno(PJ_DEFAULT_CTX);
        if (err) {
            fprintf(stderr, "Error occurred while transforming latitude=%.4f, longitude=%.4f to UTM.\n",
                    points[i].latitude, points[i].longitude);
            fprintf(stderr, "Proj error: %s\n", proj_context_errno_string(PJ_DEFAULT_CTX, err));
            return UCVM_CODE_ERROR;
        }
        point_u = xyzDest.xyzt.x;
        point_v = xyzDest.xyzt.y;

        // Point within rectangle.
        point_u -= cs241_configuration->bottom_left_corner_e;
        point_v -= cs241_configuration->bottom_left_corner_n;

        // We need to rotate that point, the number of degrees we calculated above.
        point_x = cs241_cos_rotation_angle * point_u - cs241_sin_rotation_angle * point_v;
        point_y = cs241_sin_rotation_angle * point_u + cs241_cos_rotation_angle * point_v;

        // Which point base point does that correspond to?
        load_x_coord = floor(point_x / cs241_total_width_m * (cs241_configuration->nx - 1));
        load_y_coord = floor(point_y / cs241_total_height_m * (cs241_configuration->ny - 1));

        // And on the Z-axis?
        load_z_coord = (cs241_configuration->depth / cs241_configuration->depth_interval - 1) -
                       floor(points[i].depth / cs241_configuration->depth_interval);

        // Are we outside the model's X and Y boundaries?
        if (load_x_coord > cs241_configuration->nx - 2 || load_y_coord > cs241_configuration->ny - 2 || load_x_coord < 0 || load_y_coord < 0) {
            data[i].vp = -1;
            data[i].vs = -1;
            data[i].rho = -1;
            data[i].qp = -1;
            data[i].qs = -1;
            continue;
        }

        // Get the X, Y, and Z percentages for the bilinear or cs241_trilinear interpolation below.
            double x_interval=(cs241_configuration->nx > 1) ?
                     cs241_total_width_m / (cs241_configuration->nx-1):cs241_total_width_m;
                double y_interval=(cs241_configuration->ny > 1) ?
                     cs241_total_height_m / (cs241_configuration->ny-1):cs241_total_height_m;

        x_percent = fmod(point_x, x_interval) / (x_interval);
        y_percent = fmod(point_y, y_interval) / (y_interval);
        z_percent = fmod(points[i].depth, cs241_configuration->depth_interval) / cs241_configuration->depth_interval;

        if (load_z_coord < 1) {
            // We're below the model boundaries. Bilinearly interpolate the bottom plane and use that value.
            data[i].vp = -1;
            data[i].vs = -1;
            data[i].rho = -1;
            data[i].qp = -1;
            data[i].qs = -1;

                continue;
                } else {
if ((points[i].depth < cs241_configuration->depth_interval) && (cs241_configuration->gtl == 1)) {
                           cs241_get_vs30_based_gtl(&(points[i]), &(data[i]));
                           data[i].rho=cs241_calculate_density(data[i].vs);

                        } else {
                // Read all the surrounding point properties.
                cs241_read_properties(load_x_coord,     load_y_coord,     load_z_coord,     &(surrounding_points[0]));    // Orgin.
                cs241_read_properties(load_x_coord + 1, load_y_coord,     load_z_coord,     &(surrounding_points[1]));    // Orgin + 1x
                cs241_read_properties(load_x_coord,     load_y_coord + 1, load_z_coord,     &(surrounding_points[2]));    // Orgin + 1y
                cs241_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord,     &(surrounding_points[3]));    // Orgin + x + y, forms top plane.
                cs241_read_properties(load_x_coord,     load_y_coord,     load_z_coord - 1, &(surrounding_points[4]));    // Bottom plane origin
                cs241_read_properties(load_x_coord + 1, load_y_coord,     load_z_coord - 1, &(surrounding_points[5]));    // +1x
                cs241_read_properties(load_x_coord,     load_y_coord + 1, load_z_coord - 1, &(surrounding_points[6]));    // +1y
                cs241_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord - 1, &(surrounding_points[7]));    // +x +y, forms bottom plane.

                cs241_trilinear_interpolation(x_percent, y_percent, z_percent, surrounding_points, &(data[i]));
        }
}

        // Calculate Qp and Qs.
        if (data[i].vs < 1500)
            data[i].qs = data[i].vs * 0.02;
        else
            data[i].qs = data[i].vs * 0.10;

        data[i].qp = data[i].qs * 1.5;
    }

    return SUCCESS;
}

/**
 * Retrieves the material properties (whatever is available) for the given data point, expressed
 * in x, y, and z co-ordinates.
 *
 * @param x The x coordinate of the data point.
 * @param y The y coordinate of the data point.
 * @param z The z coordinate of the data point.
 * @param data The properties struct to which the material properties will be written.
 */
void cs241_read_properties(int x, int y, int z, cs241_properties_t *data) {
    // Set everything to -1 to indicate not found.
    data->vp = -1;
    data->vs = -1;
    data->rho = -1;
    data->qp = -1;
    data->qs = -1;
    float *ptr = NULL;
    FILE *fp = NULL;
    int location = z * cs241_configuration->nx * cs241_configuration->ny + y * cs241_configuration->nx + x;

    // Check our loaded components of the model.
    if (cs241_velocity_model->vs_status == 2) {
        // Read from memory.
        ptr = (float *)cs241_velocity_model->vs;
        data->vs = ptr[location];
    } else if (cs241_velocity_model->vs_status == 1) {
        // Read from file.
        fp = (FILE *)cs241_velocity_model->vs;
        fseek(fp, location * sizeof(float), SEEK_SET);
        fread(&(data->vs), sizeof(float), 1, fp);
    }

    // Check our loaded components of the model.
    if (cs241_velocity_model->vp_status == 2) {
        // Read from memory.
        ptr = (float *)cs241_velocity_model->vp;
        data->vp = ptr[location];
    } else if (cs241_velocity_model->vp_status == 1) {
        // Read from file.
        fseek(fp, location * sizeof(float), SEEK_SET);
        fread(&(data->vp), sizeof(float), 1, fp);
    }

    // Check our loaded components of the model.
    if (cs241_velocity_model->rho_status == 2) {
        // Read from memory.
        ptr = (float *)cs241_velocity_model->rho;
        data->rho = ptr[location];
    } else if (cs241_velocity_model->rho_status == 1) {
        // Read from file.
        fseek(fp, location * sizeof(float), SEEK_SET);
        fread(&(data->rho), sizeof(float), 1, fp);
    }
}

/**
 * Trilinearly interpolates given a x percentage, y percentage, z percentage and a cube of
 * data properties in top origin format (top plane first, bottom plane second).
 *
 * @param x_percent X percentage
 * @param y_percent Y percentage
 * @param z_percent Z percentage
 * @param eight_points Eight surrounding data properties
 * @param ret_properties Returned data properties
 */
void cs241_trilinear_interpolation(double x_percent, double y_percent, double z_percent,
                             cs241_properties_t *eight_points, cs241_properties_t *ret_properties) {
    cs241_properties_t *temp_array = calloc(2, sizeof(cs241_properties_t));
    cs241_properties_t *four_points = eight_points;

    cs241_bilinear_interpolation(x_percent, y_percent, four_points, &temp_array[0]);

    // Now advance the pointer four "cs241_properties_t" spaces.
    four_points += 4;

    // Another interpolation.
    cs241_bilinear_interpolation(x_percent, y_percent, four_points, &temp_array[1]);

    // Now linearly interpolate between the two.
    cs241_linear_interpolation(z_percent, &temp_array[0], &temp_array[1], ret_properties);

    free(temp_array);
}

/**
 * Bilinearly interpolates given a x percentage, y percentage, and a plane of data properties in
 * origin, bottom-right, top-left, top-right format.
 *
 * @param x_percent X percentage.
 * @param y_percent Y percentage.
 * @param four_points Data property plane.
 * @param ret_properties Returned data properties.
 */
void cs241_bilinear_interpolation(double x_percent, double y_percent, cs241_properties_t *four_points, cs241_properties_t *ret_properties) {
    cs241_properties_t *temp_array = calloc(2, sizeof(cs241_properties_t));
    cs241_linear_interpolation(x_percent, &four_points[0], &four_points[1], &temp_array[0]);
    cs241_linear_interpolation(x_percent, &four_points[2], &four_points[3], &temp_array[1]);
    cs241_linear_interpolation(y_percent, &temp_array[0], &temp_array[1], ret_properties);
    free(temp_array);
}

/**
 * Linearly interpolates given a percentage from x0 to x1, a data point at x0, and a data point at x1.
 *
 * @param percent Percent of the way from x0 to x1 (from 0 to 1 interval).
 * @param x0 Data point at x0.
 * @param x1 Data point at x1.
 * @param ret_properties Resulting data properties.
 */
void cs241_linear_interpolation(double percent, cs241_properties_t *x0, cs241_properties_t *x1, cs241_properties_t *ret_properties) {
    ret_properties->vp  = (1 - percent) * x0->vp  + percent * x1->vp;
    ret_properties->vs  = (1 - percent) * x0->vs  + percent * x1->vs;
    ret_properties->rho = (1 - percent) * x0->rho + percent * x1->rho;
    ret_properties->qp  = (1 - percent) * x0->qp  + percent * x1->qp;
    ret_properties->qs  = (1 - percent) * x0->qs  + percent * x1->qs;
}

/**
 * Called when the model is being discarded. Free all variables.
 *
 * @return SUCCESS
 */
int cs241_finalize() {

    proj_destroy(cs241_geo2utm);
    cs241_geo2utm = NULL;

    proj_destroy(cs241_geo2aeqd);
    cs241_geo2aeqd = NULL;

    if (cs241_velocity_model) {
        free(cs241_velocity_model);
        cs241_velocity_model=NULL;
    }
    if (cs241_configuration) {
         free(cs241_configuration);
         cs241_configuration=NULL;
    }
    if (cs241_vs30_map) {
         free(cs241_vs30_map);
         cs241_vs30_map=NULL;
    }

    return SUCCESS;
}

/**
 * Returns the version information.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int cs241_version(char *ver, int len)
{
  int verlen;
  verlen = strlen(cs241_version_string);
  if (verlen > len - 1) {
    verlen = len - 1;
  }
  memset(ver, 0, len);
  strncpy(ver, cs241_version_string, verlen);
  return 0;
}

/**
 * Returns the model config information.
 *
 * @param key Config key string to return.
 * @param sz number of config terms.
 * @return Zero
 */
int cs241_config(char **config, int *sz)
{
  int len=strlen(cs241_config_string);
  if(len > 0) {
    *config=cs241_config_string;
    *sz=cs241_config_sz;
    return SUCCESS;
  }
  return FAIL;
}


/**
 * Reads the cs241_configuration file describing the various properties of CVM-S5 and populates
 * the cs241_configuration struct. This assumes cs241_configuration has been "calloc'ed" and validates
 * that each value is not zero at the end.
 *
 * @param file The cs241_configuration file location on disk to read.
 * @param config The cs241_configuration struct to which the data should be written.
 * @return Success or failure, depending on if file was read successfully.
 */
int cs241_read_configuration(char *file, cs241_configuration_t *config) {
    FILE *fp = fopen(file, "r");
    char key[40];
    char value[80];
    char line_holder[128];

    // If our file pointer is null, an error has occurred. Return fail.
    if (fp == NULL) {
        cs241_print_error("Could not open the cs241_configuration file.");
        return FAIL;
    }

    // Read the lines in the cs241_configuration file.
    while (fgets(line_holder, sizeof(line_holder), fp) != NULL) {
        if (line_holder[0] != '#' && line_holder[0] != ' ' && line_holder[0] != '\n') {
            sscanf(line_holder, "%s = %s", key, value);

            // Which variable are we editing?
            if (strcmp(key, "utm_zone") == 0)                 config->utm_zone = atoi(value);
            if (strcmp(key, "model_dir") == 0)                sprintf(config->model_dir, "%s", value);
            if (strcmp(key, "nx") == 0)                       config->nx = atoi(value);
            if (strcmp(key, "ny") == 0)                        config->ny = atoi(value);
            if (strcmp(key, "nz") == 0)                        config->nz = atoi(value);
            if (strcmp(key, "depth") == 0)                        config->depth = atof(value);
            if (strcmp(key, "top_left_corner_e") == 0)         config->top_left_corner_e = atof(value);
            if (strcmp(key, "top_left_corner_n") == 0)        config->top_left_corner_n = atof(value);
            if (strcmp(key, "top_right_corner_e") == 0)        config->top_right_corner_e = atof(value);
            if (strcmp(key, "top_right_corner_n") == 0)        config->top_right_corner_n = atof(value);
            if (strcmp(key, "bottom_left_corner_e") == 0)    config->bottom_left_corner_e = atof(value);
            if (strcmp(key, "bottom_left_corner_n") == 0)    config->bottom_left_corner_n = atof(value);
            if (strcmp(key, "bottom_right_corner_e") == 0)    config->bottom_right_corner_e = atof(value);
            if (strcmp(key, "bottom_right_corner_n") == 0)    config->bottom_right_corner_n = atof(value);
            if (strcmp(key, "depth_interval") == 0)            config->depth_interval = atof(value);
                        if (strcmp(key, "p0") == 0)                                             config->p0 = atof(value);
                        if (strcmp(key, "p1") == 0)                                             config->p1 = atof(value);
                        if (strcmp(key, "p2") == 0)                                             config->p2 = atof(value);
                        if (strcmp(key, "p3") == 0)                                             config->p3 = atof(value);
                        if (strcmp(key, "p4") == 0)                                             config->p4 = atof(value);
                        if (strcmp(key, "p5") == 0)                                             config->p5 = atof(value);
                        if (strcmp(key, "gtl") == 0) {
                                if (strcmp(value, "on") == 0) config->gtl = 1;
                                else config->gtl = 0;
                        }
// anything else, just ignore
        }
    }

    // Have we set up all cs241_configuration parameters?
    if (config->utm_zone == 0 || config->nx == 0 || config->ny == 0 || config->nz == 0 || config->model_dir[0] == '\0' ||
        config->top_left_corner_e == 0 || config->top_left_corner_n == 0 || config->top_right_corner_e == 0 ||
        config->top_right_corner_n == 0 || config->bottom_left_corner_e == 0 || config->bottom_left_corner_n == 0 ||
        config->bottom_right_corner_e == 0 || config->bottom_right_corner_n == 0 || config->depth == 0 ||
        config->depth_interval == 0) {
        cs241_print_error("One cs241_configuration parameter not specified. Please check your cs241_configuration file.");
        return FAIL;
    }

    fclose(fp);

    return SUCCESS;
}

/**
 *  * Calculates the density based off of Vs. Based on Nafe-Drake scaling relationship.
 *   *
 *    * @param vs The Vs value off which to scale.
 *     * @return Density, in g/m^3.
 *      */
double cs241_calculate_density(double vs) {
        double retVal;
        vs = vs / 1000;
        retVal = cs241_configuration->p0 + cs241_configuration->p1 * vs + cs241_configuration->p2 * pow(vs, 2) +
                         cs241_configuration->p3 * pow(vs, 3) + cs241_configuration->p4 * pow(vs, 4) + cs241_configuration->p5 * pow(vs, 5);
        retVal = retVal * 1000;
        return retVal;
}


/**
 * Prints the error string provided.
 *
 * @param err The error string to print out to stderr.
 */
void cs241_print_error(char *err) {
    fprintf(stderr, "An error has occurred while executing CS241. The error was:\n\n");
    fprintf(stderr, "%s", err);
    fprintf(stderr, "\n\nPlease contact software@scec.org and describe both the error and a bit\n");
    fprintf(stderr, "about the computer you are running CS241 on (Linux, Mac, etc.).\n");
}

/**
 * Tries to read the model into memory.
 *
 * @param model The model parameter struct which will hold the pointers to the data either on disk or in memory.
 * @return 2 if all files are read to memory, SUCCESS if file is found but at least 1
 * is not in memory, FAIL if no file found.
 */
int cs241_try_reading_model(cs241_model_t *model) {
    double base_malloc = cs241_configuration->nx * cs241_configuration->ny * cs241_configuration->nz * sizeof(float);
    int file_count = 0;
    int all_read_to_memory = 1;
    char current_file[128];
    FILE *fp;

    // Let's see what data we actually have.
    sprintf(current_file, "%s/vp.dat", cs241_iteration_directory);
    if (access(current_file, R_OK) == 0) {
        model->vp = malloc(base_malloc);
        if (model->vp != NULL) {
            // Read the model in.
            fp = fopen(current_file, "rb");
            fread(model->vp, 1, base_malloc, fp);
            fclose(fp);
            model->vp_status = 2;
        } else {
            all_read_to_memory = 0;
            model->vp = fopen(current_file, "rb");
            model->vp_status = 1;
        }
        file_count++;
    }

    sprintf(current_file, "%s/vs.dat", cs241_iteration_directory);
    if (access(current_file, R_OK) == 0) {
        model->vs = malloc(base_malloc);
        if (model->vs != NULL) {
            // Read the model in.
            fp = fopen(current_file, "rb");
            fread(model->vs, 1, base_malloc, fp);
            fclose(fp);
            model->vs_status = 2;
        } else {
            all_read_to_memory = 0;
            model->vs = fopen(current_file, "rb");
            model->vs_status = 1;
        }
        file_count++;
    }

    sprintf(current_file, "%s/density.dat", cs241_iteration_directory);
    if (access(current_file, R_OK) == 0) {
        model->rho = malloc(base_malloc);
        if (model->rho != NULL) {
            // Read the model in.
            fp = fopen(current_file, "rb");
            fread(model->rho, 1, base_malloc, fp);
            fclose(fp);
            model->rho_status = 2;
        } else {
            all_read_to_memory = 0;
            model->rho = fopen(current_file, "rb");
            model->rho_status = 1;
        }
        file_count++;
    }

    sprintf(current_file, "%s/qp.dat", cs241_iteration_directory);
    if (access(current_file, R_OK) == 0) {
        model->qp = malloc(base_malloc);
        if (model->qp != NULL) {
            // Read the model in.
            fp = fopen(current_file, "rb");
            fread(model->qp, 1, base_malloc, fp);
            fclose(fp);
            model->qp_status = 2;
        } else {
            all_read_to_memory = 0;
            model->qp = fopen(current_file, "rb");
            model->qp_status = 1;
        }
        file_count++;
    }

    sprintf(current_file, "%s/qs.dat", cs241_iteration_directory);
    if (access(current_file, R_OK) == 0) {
        model->qs = malloc(base_malloc);
        if (model->qs != NULL) {
            // Read the model in.
            fp = fopen(current_file, "rb");
            fread(model->qs, 1, base_malloc, fp);
            fclose(fp);
            model->qs_status = 2;
        } else {
            all_read_to_memory = 0;
            model->qs = fopen(current_file, "rb");
            model->qs_status = 1;
        }
        file_count++;
    }

    if (file_count == 0)
        return FAIL;
    else if (file_count > 0 && all_read_to_memory == 0)
        return SUCCESS;
    else
        return 2;
}


/**
 * Reads the format of the Vs30 data e-tree. This file location is typically specified
 * in the cs241_configuration file of the model.
 *
 * @param filename The e-tree's file location from which to read.
 * @param map The outputted map cs241_configuration structure.
 */
int cs241_read_vs30_map(char *filename, cs241_vs30_map_config_t *map) {
    char appmeta[512];
    char *token;
    int index = 0, retVal = 0;
    map->vs30_map = etree_open(filename, O_RDONLY, 64, 0, 3);
    retVal = snprintf(appmeta, sizeof(appmeta), "%s", etree_getappmeta(map->vs30_map));

    if (retVal >= 0 && retVal < 128) {
        return FAIL;
    }

    // Now we need to parse the map cs241_configuration.
    index = 0;
    token = strtok(appmeta, "|");

    while (token != NULL) {
        switch (index) {
        case 0:
            snprintf(map->type, sizeof(map->type), "%s", token);
            break;
        case 1:
            snprintf(map->description, sizeof(map->description), "%s", token);
            break;
        case 2:
            snprintf(map->author, sizeof(map->author), "%s", token);
            break;
        case 3:
            snprintf(map->date, sizeof(map->date), "%s", token);
            break;
        case 4:
            sscanf(token, "%lf", &(map->spacing));
            break;
        case 5:
            snprintf(map->schema, sizeof(map->schema), "%s", token);
            break;
        case 6:
            snprintf(map->projection, sizeof(map->projection), "%s", token);
            break;
        case 7:
            sscanf(token, "%lf,%lf,%lf", &(map->origin_point.longitude), &(map->origin_point.latitude),
                    &(map->origin_point.depth));
            break;
        case 8:
            sscanf(token, "%lf", &(map->rotation));
            break;
        case 9:
            sscanf(token, "%lf,%lf,%lf", &(map->x_dimension), &(map->y_dimension), &(map->z_dimension));
            break;
        case 10:
            sscanf(token, "%u,%u,%u", &(map->x_ticks), &(map->y_ticks), &(map->z_ticks));
            break;
        default:
            fprintf(stderr, "Unexpected metadata. Please check your Vs30 e-tree within UCVM.\n");
            return FAIL;
            break;
        }
        index++;
        token = strtok(NULL, "|");
    }

    return SUCCESS;

}

/**
 * Given a latitude and longitude in WGS84 co-ordinates, we find the corresponding e-tree octant
 * in the Vs30 map e-tree and read the value as well as interpolate bilinearly.
 *
 * @param longitude The longitude in WGS84 format.
 * @param latitude The latitude in WGS84 format.
 * @param map The Vs30 map structure as defined during the initialization procedure.
 * @return The Vs30 value at that point, or -1 if outside the boundaries.
 */
double cs241_get_vs30_value(double longitude, double latitude, cs241_vs30_map_config_t *map) {
    // Convert both points to UTM.

    double point_x, point_y;
    double origin_x, origin_y;

    double temp_rotated_point_x = 0.0, temp_rotated_point_y = 0.0;
    double rotated_point_x = 0.0, rotated_point_y = 0.0;

    double percent = 0.0;

    int loc_x = 0, loc_y = 0;
    etree_addr_t addr;
    cs241_vs30_mpayload_t vs30_payload[4];

    int max_level = ceil(log(map->x_dimension / map->spacing) / log(2.0));

    etree_tick_t edgetics = (etree_tick_t)1 << (ETREE_MAXLEVEL - max_level);
    double map_edgesize = map->x_dimension / (double)((etree_tick_t)1<<max_level);

    PJ_COORD xyzSrc = proj_coord(latitude, longitude, 0.0, HUGE_VAL);
    PJ_COORD xyzDest = proj_trans(cs241_geo2aeqd, PJ_FWD, xyzSrc);
    point_x = xyzDest.xyzt.x;
    point_y = xyzDest.xyzt.y;
    fprintf(stderr,"  xyzDest (%.6f)  (%.6f)\n",point_x, point_y);

    xyzSrc = proj_coord(map->origin_point.latitude, map->origin_point.longitude, 0.0, HUGE_VAL);
    xyzDest = proj_trans(cs241_geo2aeqd, PJ_FWD, xyzSrc);
    origin_x = xyzDest.xyzt.x;
    origin_y = xyzDest.xyzt.y;
    fprintf(stderr,"  vs30Dest (%.6f)  (%.6f)\n",origin_x, origin_y);

    // Now that both are in UTM, we can subtract and rotate.
    temp_rotated_point_x = point_x - origin_x;
    temp_rotated_point_y = point_y - origin_y;

    rotated_point_x = cs241_cos_vs30_rotation_angle * temp_rotated_point_x - cs241_sin_vs30_rotation_angle * temp_rotated_point_y;
    rotated_point_y = cs241_sin_vs30_rotation_angle * temp_rotated_point_x + cs241_cos_vs30_rotation_angle * temp_rotated_point_y;

    // Are we within the box?
    if (rotated_point_x < 0 || rotated_point_y < 0 || rotated_point_x > map->x_dimension ||
        rotated_point_y > map->y_dimension) return -1;

    // Get the integer location of the grid point within the map.
    loc_x = floor(rotated_point_x / map_edgesize);
    loc_y = floor(rotated_point_y / map_edgesize);

    // We need the four surrounding points for bilinear interpolation.
    addr.level = ETREE_MAXLEVEL;
    addr.x = loc_x * edgetics; addr.y = loc_y * edgetics; addr.z = 0;
    /* Adjust addresses for edges of grid */
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
    etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[0]));
    addr.x = (loc_x + 1) * edgetics; addr.y = loc_y * edgetics;
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
    etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[1]));
    addr.x = loc_x * edgetics; addr.y = (loc_y + 1) * edgetics;
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
    etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[2]));
    addr.x = (loc_x + 1) * edgetics; addr.y = (loc_y + 1) * edgetics;
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
    etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[3]));

    percent = fmod(rotated_point_x / map->spacing, map->spacing) / map->spacing;
    vs30_payload[0].vs30 = percent * vs30_payload[0].vs30 + (1 - percent) * vs30_payload[1].vs30;
    vs30_payload[1].vs30 = percent * vs30_payload[2].vs30 + (1 - percent) * vs30_payload[3].vs30;

    return vs30_payload[0].vs30;
}

/**
 * Gets the GTL value using the Wills and Wald dataset, given a latitude, longitude and depth.
 *
 * @param point The point at which to retrieve the property. Note, depth is ignored.
 * @param data The material properties at the point specified, or -1 if not found.
 * @return Success or failure.
 */
int cs241_get_vs30_based_gtl(cs241_point_t *point, cs241_properties_t *data) {
        double a = 0.5, b = 0.6, c = 0.5;
    double percent_z = point->depth / cs241_configuration->depth_interval;
    double f = 0.0, g = 0.0;
    double vs30 = 0.0, vp30 = 0.0;

    // Double check that we're above the first layer.
    if (percent_z > 1) return FAIL;

    // Query for the point at depth_interval.
    cs241_point_t *pt = calloc(1, sizeof(cs241_point_t));
    cs241_properties_t *dt = calloc(1, sizeof(cs241_properties_t));

    pt->latitude = point->latitude;
    pt->longitude = point->longitude;
    pt->depth = cs241_configuration->depth_interval;

    if (cs241_query(pt, dt, 1) != SUCCESS) return FAIL;

    // Now we need the Vs30 data value.
    vs30 = cs241_get_vs30_value(point->longitude, point->latitude, cs241_vs30_map);

    if (vs30 == -1) {
        data->vp = -1;
        data->vs = -1;
    } else {
        // Get the point's material properties within the GTL.
        f = percent_z + b * (percent_z - pow(percent_z, 2.0f));
        g = a - a * percent_z + c * (pow(percent_z, 2.0f) + 2.0 * sqrt(percent_z) - 3.0 * percent_z);
        data->vs = f * dt->vs + g * vs30;
//fprintf(stderr,"XXX f %f and g %f\n", f, g);
        vs30 = vs30 / 1000;
        vp30 = 0.9409 + 2.0947 * vs30 - 0.8206 * pow(vs30, 2.0f) + 0.2683 * pow(vs30, 3.0f) - 0.0251 * pow(vs30, 4.0f);
        vp30 = vp30 * 1000;
        data->vp = f * dt->vp + g * vp30;
    }

    free(pt);
    free(dt);

    return SUCCESS;
}


// The following functions are for dynamic library mode. If we are compiling
// a static library, these functions must be disabled to avoid conflicts.
#ifdef DYNAMIC_LIBRARY

/**
 * Init function loaded and called by the UCVM library. Calls cs241_init.
 *
 * @param dir The directory in which UCVM is installed.
 * @return Success or failure.
 */
int model_init(const char *dir, const char *label) {
    return cs241_init(dir, label);
}

/**
 * Query function loaded and called by the UCVM library. Calls cs241_query.
 *
 * @param points The basic_point_t array containing the points.
 * @param data The basic_properties_t array containing the material properties returned.
 * @param numpoints The number of points in the array.
 * @return Success or fail.
 */
int model_query(cs241_point_t *points, cs241_properties_t *data, int numpoints) {
    return cs241_query(points, data, numpoints);
}

/**
 * Finalize function loaded and called by the UCVM library. Calls cs241_finalize.
 *
 * @return Success
 */
int model_finalize() {
    return cs241_finalize();
}

/**
 * Version function loaded and called by the UCVM library. Calls cs241_version.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int model_version(char *ver, int len) {
        return cs241_version(ver, len);
}

/**
 * Version function loaded and called by the UCVM library. Calls cs241_config.
 *
 * @param config Config string to return.
 * @param sz number of config terms
 * @return Zero
 */
int model_config(char **config, int *sz) {
    return cs241_config(config, sz);
}

int (*get_model_init())(const char *, const char *) {
    return &cs241_init;
}
int (*get_model_query())(cs241_point_t *, cs241_properties_t *, int) {
    return &cs241_query;
}
int (*get_model_finalize())() {
    return &cs241_finalize;
}
int (*get_model_version())(char *, int) {
    return &cs241_version;
}
int (*get_model_config())(char **, int*) {
    return &cs241_config;
}

#endif