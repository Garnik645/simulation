#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <stdbool.h>
#include <math.h>

#include "cJSON.h"

#ifdef _WIN32
#define PATH_SEPARATOR '\\'
#else
#define PATH_SEPARATOR '/'
#endif

typedef long double ld;

const int thread_num = 12;
const ld derivative_delta = 0.0001;

char *out_dir;

int n_steps;
int period;

ld delta_time;
ld cutoff;
ld a;
ld b;
ld c;

typedef struct ParticleType {
  char *name;
  ld m;
  ld sigma;
  ld eps;
} ParticleType;
int particle_types_size;
ParticleType *particle_types;

typedef struct Particle {
  ParticleType type;
  char *type_name;
  ld x;
  ld y;
  ld z;
  ld vx;
  ld vy;
  ld vz;
  int index;
} Particle;
int particles_size;
Particle *particles;

ld **eps;
ld **sigma;

typedef struct Vector {
  ld x;
  ld y;
  ld z;
} Vector;

// CMD -i /path/to/input.csv -p /path/to/params.json -n NSTEPS -period PERIOD -o TRAJ_DIR

char *get_file_content(const char *path) {
  FILE *file = fopen(path, "r");
  if (file == NULL) {
    perror("Failed to open file");
    exit(EXIT_FAILURE);
  }

  // Determine the file size
  fseek(file, 0, SEEK_END);
  long file_size = ftell(file);
  rewind(file);

  // Allocate memory for the file contents
  char *content = (char *) malloc(file_size + 1);
  if (content == NULL) {
    perror("Failed to allocate memory");
    fclose(file);
    exit(EXIT_FAILURE);
  }

  // Read the file contents into the buffer
  size_t bytes_read = fread(content, 1, file_size, file);
  content[bytes_read] = '\0'; // Null-terminate the string

  // Close the file
  fclose(file);
  return content;
}

void get_csv_data(const char *csv_path) {
  char *content = get_file_content(csv_path);
  char *result = strchr(content, '\n');
  int number_of_lines = 1;
  for (char *curr = result; *curr != '\0'; ++curr) {
    if (*curr == ',') {
      *curr = ' ';
    } else if (*curr == '\n') {
      ++number_of_lines;
    }
  }
  particles = (Particle *) malloc(sizeof(Particle) * number_of_lines);
  int i = 0;
  while (*result != '\0') {
    particles[i].type_name = malloc(sizeof(char) * 10);
    int ret = sscanf(result, "%d %s %Lf %Lf %Lf %Lf %Lf %Lf", &(particles[i].index), particles[i].type_name,
                     &(particles[i].x),
                     &(particles[i].y),
                     &(particles[i].z), &(particles[i].vx), &(particles[i].vy), &(particles[i].vz));
    if (ret != 8) {
      break;
    }
    ++i;
    ++result;
    while (*result != '\n' && *result != '\0') {
      ++result;
    }
  }
  free(content);
  particles_size = i;
  for (i = 0; i < particles_size; ++i) {
    bool found = false;
    for (int j = 0; j < particle_types_size; ++j) {
      if (strcmp(particles[i].type_name, particle_types[j].name) == 0) {
        particles[i].type = particle_types[j];
        found = true;
        break;
      }
    }
    if (!found) {
      printf("Particle type not found!");
      exit(EXIT_FAILURE);
    }
  }
}

void parse_args(int argc, char *argv[]) {
  assert(argc == 11);
  assert(strcmp(argv[1], "-i") == 0);
  assert(strcmp(argv[3], "-p") == 0);
  assert(strcmp(argv[5], "-n") == 0);
  assert(strcmp(argv[7], "-period") == 0);
  assert(strcmp(argv[9], "-o") == 0);
  const char *input_csv_path = argv[2];
  const char *params_json_path = argv[4];
  n_steps = atoi(argv[6]);
  period = atoi(argv[8]);
  out_dir = argv[10];

  char *params_content = get_file_content(params_json_path);
  cJSON *params = cJSON_Parse(params_content);
  free(params_content);
  if (params == NULL) {
    const char *error_ptr = cJSON_GetErrorPtr();
    if (error_ptr != NULL) {
      fprintf(stderr, "Error before: %s\n", error_ptr);
    }
    cJSON_Delete(params);
    exit(EXIT_FAILURE);
  }

  cJSON *dt_json = cJSON_GetObjectItem(params, "dt");
  if (dt_json == NULL || !cJSON_IsNumber(dt_json)) {
    exit(EXIT_FAILURE);
  }
  delta_time = dt_json->valuedouble;

  cJSON *cutoff_json = cJSON_GetObjectItem(params, "cutoff");
  if (cutoff_json == NULL || !cJSON_IsNumber(cutoff_json)) {
    exit(EXIT_FAILURE);
  }
  cutoff = cutoff_json->valuedouble;

  cJSON *pbc_json = cJSON_GetObjectItem(params, "pbc");
  if (pbc_json == NULL || !cJSON_IsArray(pbc_json)) {
    exit(EXIT_FAILURE);
  }
  cJSON *element;
  ld box[3];
  int i = 0;
  cJSON_ArrayForEach(element, pbc_json) {
    // Access each object in the array
    box[i] = element->valuedouble;
    ++i;
  }
  a = box[0];
  b = box[1];
  c = box[2];

  cJSON *particle_types_json = cJSON_GetObjectItem(params, "particle_types");
  if (particle_types_json == NULL && !cJSON_IsObject(particle_types_json)) {
    exit(EXIT_FAILURE);
  }
  particle_types_size = cJSON_GetArraySize(particle_types_json);
  particle_types = (ParticleType *) malloc(sizeof(ParticleType) * particle_types_size);
  i = 0;
  cJSON_ArrayForEach(element, particle_types_json) {
    particle_types[i].name = strdup(element->string);
    particle_types[i].m = cJSON_GetObjectItem(element, "m")->valuedouble;
    particle_types[i].sigma = cJSON_GetObjectItem(element, "lj_sigma")->valuedouble;
    particle_types[i].eps = cJSON_GetObjectItem(element, "lj_eps")->valuedouble;
    ++i;
  }
  cJSON_Delete(params);

  get_csv_data(input_csv_path);
}

void calc_eps_sigma() {
  eps = (ld **) malloc(sizeof(ld *) * particles_size);
  sigma = (ld **) malloc(sizeof(ld *) * particles_size);
  for (int i = 0; i < particles_size; ++i) {
    eps[i] = (ld *) malloc(sizeof(ld) * particles_size);
    sigma[i] = (ld *) malloc(sizeof(ld) * particles_size);
  }
  for (int i = 0; i < particles_size; ++i) {
    for (int j = 0; j < particles_size; ++j) {
      if (i != j) {
        eps[i][j] = sqrt(particles[i].type.eps * particles[j].type.eps);
        sigma[i][j] = (particles[i].type.sigma + particles[j].type.sigma) / 2.0;
      }
    }
  }
}

ld min(ld x, ld y) {
  return (x < y) ? x : y;
}

ld Vlj(int from, int to, ld tox, ld toy, ld toz) {
  ld dx = particles[from].x - tox;
  dx = min(dx, a - dx);
  ld dy = particles[from].y - toy;
  dy = min(dy, b - dy);
  ld dz = particles[from].z - toz;
  dz = min(dz, c - dz);
  ld dist2 = dx * dx + dy * dy + dz * dz;
  if (dist2 > cutoff * cutoff) {
    return 0.0;
  }
  ld tmp = pow(sigma[from][to], 6) / pow(dist2, 3);
  ld result = 4 * eps[from][to] * tmp * (tmp - 1);
  return result;
}

ld get_derivative_x(int from, int to) {
  ld lhs = Vlj(from, to, particles[to].x + derivative_delta, particles[to].y, particles[to].z);
  ld rhs = Vlj(from, to, particles[to].x - derivative_delta, particles[to].y, particles[to].z);
  return (lhs - rhs) / (2 * derivative_delta);
}

ld get_derivative_y(int from, int to) {
  ld lhs = Vlj(from, to, particles[to].x, particles[to].y + derivative_delta, particles[to].z);
  ld rhs = Vlj(from, to, particles[to].x, particles[to].y - derivative_delta, particles[to].z);
  return (lhs - rhs) / (2 * derivative_delta);
}

ld get_derivative_z(int from, int to) {
  ld lhs = Vlj(from, to, particles[to].x, particles[to].y, particles[to].z + derivative_delta);
  ld rhs = Vlj(from, to, particles[to].x, particles[to].y, particles[to].z - derivative_delta);
  return (lhs - rhs) / (2 * derivative_delta);
}

Vector affect(int from, int to) {
  Vector result;
  result.x = -get_derivative_x(from, to);
  result.y = -get_derivative_y(from, to);
  result.z = -get_derivative_z(from, to);
  return result;
}

typedef struct args {
  Vector *force;
  int index;
} args;

void *get_single_force(void *arg) {
  args *arg_t = (args *)arg;
  Vector *force = arg_t->force;
  int index = arg_t->index;
  for (int i = index; i < particles_size; i += thread_num) {
    for (int j = 0; j < particles_size; ++j) {
      if (i != j) {
        Vector new_force = affect(j, i);
        force[i].x += new_force.x;
        force[i].y += new_force.y;
        force[i].z += new_force.z;
      }
    }
  }
  return NULL;
}

Vector *get_force() {
  Vector *force = (Vector *) malloc(sizeof(Vector) * particles_size);
  for (int i = 0; i < particles_size; ++i) {
    force[i].x = 0;
    force[i].y = 0;
    force[i].z = 0;
  }
  pthread_t* threads = (pthread_t *)malloc(sizeof(pthread_t) * thread_num);
  args *thread_args = (args *)malloc(sizeof(args) * thread_num);
  for(int i = 0; i < thread_num; ++i)
  {
    thread_args[i].force = force;
    thread_args[i].index = i;
    int thread_out = pthread_create(&threads[i], NULL, get_single_force, &thread_args[i]);
    // check for errors
    if(thread_out != 0)
    {
      printf("Error while creating thread!\n");
      exit(thread_out);
    }
  }
  for(int i = 0; i < thread_num; ++i)
  {
    void *thread_result;
    // join thread
    int thread_out = pthread_join(threads[i], &thread_result);

    // check for errors
    if(thread_out != 0)
    {
      printf("Error while joining thread!\n");
      exit(thread_out);
    }
  }
  free(thread_args);
  free(threads);
  return force;
}

void move(Vector *force) {
  for (int i = 0; i < particles_size; ++i) {
    Vector acceleration;
    acceleration.x = force[i].x / particles[i].type.m;
    acceleration.y = force[i].y / particles[i].type.m;
    acceleration.z = force[i].z / particles[i].type.m;

    particles[i].vx += acceleration.x * delta_time;
    particles[i].vy += acceleration.y * delta_time;
    particles[i].vz += acceleration.z * delta_time;

    particles[i].x += particles[i].vx * delta_time + acceleration.x * delta_time * delta_time / 2;
    while (particles[i].x < 0) {
      particles[i].x += a;
    }
    while (particles[i].x > a) {
      particles[i].x -= a;
    }

    particles[i].y += particles[i].vy * delta_time + acceleration.y * delta_time * delta_time / 2;
    while (particles[i].y < 0) {
      particles[i].y += b;
    }
    while (particles[i].y > b) {
      particles[i].y -= b;
    }

    particles[i].z += particles[i].vz * delta_time + acceleration.z * delta_time * delta_time / 2;
    while (particles[i].z < 0) {
      particles[i].z += c;
    }
    while (particles[i].z > c) {
      particles[i].z -= c;
    }
  }
}

void join_path(char *result, const char *dir, const char *filename) {
  snprintf(result, 256, "%s%c%s", dir, PATH_SEPARATOR, filename);
}

void get_file_path(char *result, int index) {
  char filename[256];
  char numbers[10];
  for (int i = 5; i >= 0; --i) {
    numbers[i] = (index % 10) + '0';
    index /= 10;
  }
  numbers[6] = '\0';
  char *file_prefix = "step_";
  char *file_suffix = ".csv";
  *filename = '\0';
  strcat(filename, file_prefix);
  strcat(filename, numbers);
  strcat(filename, file_suffix);
  join_path(result, out_dir, filename);
}

void write_frame(int index, Vector *force) {
  char file_path[256];
  get_file_path(file_path, index);

  FILE *file = fopen(file_path, "w");
  if (file == NULL) {
    printf("Error opening file.\n");
    exit(EXIT_FAILURE);
  }

  fprintf(file, "particle_type,x,y,z,vx,vy,vz,fx,fy,fz\n");
  for (int i = 0; i < particles_size; ++i) {
    fprintf(file, "%s,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf\n", particles[i].type_name, particles[i].x, particles[i].y,
            particles[i].z,
            particles[i].vx, particles[i].vy, particles[i].vz, force[i].x, force[i].y, force[i].z);
  }
  fclose(file);
}

void iterate() {
  for (int i = 0; i <= n_steps; ++i) {
    Vector *force = get_force();
    move(force);
    if (i % period == 0) {
      write_frame(i, force);
    }
    free(force);
  }
}

int main(int argc, char *argv[]) {
  parse_args(argc, argv);
  calc_eps_sigma();
  iterate();
}