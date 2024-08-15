/**
 * @file cs248.h
 * @brief Main header file for cs248 library.
 * @version 1.0
 *
 * Delivers the cs248 model 
 * base on original cs248
 *
 */

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "proj.h"

// Constants
#ifndef M_PI
	/** Defines pi */
	#define M_PI 3.14159265358979323846
#endif
#define DEG_TO_RAD M_PI / 180.0

/** Defines a return value of success */
#define SUCCESS 0
/** Defines a return value of failure */
#define FAIL 1

/* config string */
#define CS248_CONFIG_MAX 1000

// Structures
/** Defines a point (latitude, longitude, and depth) in WGS84 format */
typedef struct cs248_point_t {
	/** Longitude member of the point */
	double longitude;
	/** Latitude member of the point */
	double latitude;
	/** Depth member of the point */
	double depth;
} cs248_point_t;

/** Defines the material properties this model will retrieve. */
typedef struct cs248_properties_t {
	/** P-wave velocity in meters per second */
	double vp;
	/** S-wave velocity in meters per second */
	double vs;
	/** Density in g/m^3 */
	double rho;
	/** Qp */
	double qp;
	/** Qs */
	double qs;
} cs248_properties_t;

/** The CVM-S5 configuration structure. */
typedef struct cs248_configuration_t {
	/** The zone of UTM projection */
	int utm_zone;
	/** The model directory */
	char model_dir[128];
	/** Number of x points */
	int nx;
	/** Number of y points */
	int ny;
	/** Number of z points */
	int nz;
	/** Depth in meters */
	double depth;
	/** Top left corner easting in UTM projection */
	double top_left_corner_e;
	/** Top left corner northing in UTM projection */
	double top_left_corner_n;
	/** Top right corner easting in UTM projection */
	double top_right_corner_e;
	/** Top right corner northing in UTM projection */
	double top_right_corner_n;
	/** Bottom left corner easting in UTM projection */
	double bottom_left_corner_e;
	/** Bottom left corner northing in UTM projection */
	double bottom_left_corner_n;
	/** Bottom right corner easting in UTM projection */
	double bottom_right_corner_e;
	/** Bottom right corner northing in UTM projection */
	double bottom_right_corner_n;
	/** Z interval for the data */
	double depth_interval;
        /** The data access seek method, fast-X, or fast-Y */
        char seek_axis[128];
        /** The data seek direction, bottom-up, or top-down */
        char seek_direction[128];
        /** trilinear interploation; */
        int interpolation;
} cs248_configuration_t;

/** The model structure which points to available portions of the model. */
typedef struct cs248_model_t {
	/** A pointer to the Vs data either in memory or disk. Null if does not exist. */
	void *vs;
	/** Vs status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int vs_status;
	/** A pointer to the Vp data either in memory or disk. Null if does not exist. */
	void *vp;
	/** Vp status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int vp_status;
	/** A pointer to the rho data either in memory or disk. Null if does not exist. */
	void *rho;
	/** Rho status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int rho_status;
	/** A pointer to the Qp data either in memory or disk. Null if does not exist. */
	void *qp;
	/** Qp status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int qp_status;
	/** A pointer to the Qs data either in memory or disk. Null if does not exist. */
	void *qs;
	/** Qs status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int qs_status;
} cs248_model_t;

// UCVM API Required Functions
#ifdef DYNAMIC_LIBRARY

/** Initializes the model */
int model_init(const char *dir, const char *label);
/** Cleans up the model (frees memory, etc.) */
int model_finalize();
/** Returns version information */
int model_version(char *ver, int len);
/** Queries the model */
int model_query(cs248_point_t *points, cs248_properties_t *data, int numpts);
int model_config(char **config, int *sz);

int (*get_model_init())(const char *, const char *);
int (*get_model_query())(cs248_point_t *, cs248_properties_t *, int);
int (*get_model_finalize())();
int (*get_model_version())(char *, int);
int (*get_model_config())(char **, int*);

#endif

// cs248 Related Functions

/** Initializes the model */
int cs248_init(const char *dir, const char *label);
/** Cleans up the model (frees memory, etc.) */
int cs248_finalize();
/** Returns version information */
int cs248_version(char *ver, int len);
/** Queries the model */
int cs248_query(cs248_point_t *points, cs248_properties_t *data, int numpts);

// Non-UCVM Helper Functions
/** Reads the configuration file. */
int cs248_read_configuration(char *file, cs248_configuration_t *config);
int cs248_dump_configuration(cs248_configuration_t *config);
/** Prints out the error string. */
void cs248_print_error(char *err);
/** Retrieves the value at a specified grid point in the model. */
void cs248_read_properties(int x, int y, int z, cs248_properties_t *data);
/** Attempts to malloc the model size in memory and read it in. */
int cs248_try_reading_model(cs248_model_t *model);

// Interpolation Functions
/** Linearly interpolates two cs248_properties_t structures */
void cs248_linear_interpolation(double percent, cs248_properties_t *x0, cs248_properties_t *x1, cs248_properties_t *ret_properties);
/** Bilinearly interpolates the properties. */
void cs248_bilinear_interpolation(double x_percent, double y_percent, cs248_properties_t *four_points, cs248_properties_t *ret_properties);
/** Trilinearly interpolates the properties. */
void cs248_trilinear_interpolation(double x_percent, double y_percent, double z_percent, cs248_properties_t *eight_points,
							 cs248_properties_t *ret_properties);
