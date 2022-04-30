#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <pthread.h>

const int BIN_COUNT = 256;
struct args {
    int start;
    int length;
    int* histogram;
    int* image;
};

void die(const char* msg){
    if (errno != 0)
        perror(msg);
    else
        fprintf(stderr, "error: %s\n", msg);
    exit(1);
}

void* update_histogram(void* thread_args){
    struct args* args = (struct args*)thread_args;
    int start = args->start;
    int end = start + args->length;
    int* image = args->image;
    int* histogram = args->histogram;
    int* local_histogram = (int*)calloc(BIN_COUNT, sizeof(int));

    for (int i = start; i < end; ++i)
        local_histogram[image[i]] += 1;

    for (int i = 0; i < BIN_COUNT; ++i) {
        if (!local_histogram[i])
            continue;
        __sync_fetch_and_add(&histogram[i], local_histogram[i]);
    }
    return NULL;
}

void generate_image(int num_rows, int num_cols, int* image){
    for (int i = 0; i < num_cols * num_rows; ++i)
        image[i] = rand() % 256; // 255 + 1 for num bins
}

void read_image(const char* image_path, int num_rows, int num_cols, int* image){
    char format[3];
    FILE* f;
    unsigned imgw, imgh, maxv, v;
    size_t i;

    printf("Reading PGM data from %s...\n", image_path);

    if (!(f = fopen(image_path, "r")))
        die("fopen");

    fscanf(f, "%2s", format);
    if (format[0] != 'P' || format[1] != '2')
        die("only ASCII PGM input is supported");

    if (fscanf(f, "%u", &imgw) != 1 ||
        fscanf(f, "%u", &imgh) != 1 ||
        fscanf(f, "%u", &maxv) != 1)
        die("invalid input");

    if (imgw != num_cols || imgh != num_rows){
        fprintf(stderr, "input data size (%ux%u) does not match cylinder size (%zux%zu)\n",
            imgw, imgh, num_cols, num_rows);
        die("invalid input");
    }

    for (i = 0; i < num_cols * num_rows; ++i){
        if (fscanf(f, "%u", &v) != 1)
            die("invalid data");
        image[i] = ((int)v * 255) / maxv; // 255 for num bins
    }
    fclose(f);
}

void print_histo(int* histo){
    for (int i = 0; i < 256; ++i){
        if (i != 0 && (i % 10 == 0))
            printf("\n");
        printf("%d ", histo[i]);
    }
    printf("\n");
}

void print_image(int num_rows, int num_cols, int* image){
    int index = 0;
    for (int i = 0; i < num_rows; ++i){
        for (int j = 0; j < num_cols; ++j){
            index = i * num_cols + j;
            printf("%d ", image[index]);
        }
    }
    printf("\n");
}

void histogram(int* histo, int* image, int num_threads, int image_size){
    pthread_t* threads = (pthread_t*)malloc(sizeof(pthread_t) * num_threads);

    int task_size = image_size / num_threads;
    int last_task_size = image_size - (num_threads - 1) * task_size;
    int task_start = 0;

    for (int i = 0; i < num_threads; ++i){
        struct args* args = (struct args*)malloc(sizeof(struct args));
        args->start = task_start;
        args->length = i == num_threads - 1 ? last_task_size : task_size;
        args->histogram = histo;
        args->image = image;

        pthread_create(&threads[i], NULL, &update_histogram, args);
        task_start += task_size;
    }

    for (int i = 0; i < num_threads; ++i){
        pthread_join(threads[i], NULL);
    }
}

int main(int argc, char* argv[]){
    int c;
    int seed = 42;
    const char* image_path = 0;
    image_path = "./input/big_2000x2000.pgm";
    int gen_image = 0;
    int debug = 0;

    int num_rows = 2000;
    int num_cols = 2000;
    int num_threads = 1;

    struct timespec before, after;

    int* histo = (int*)calloc(256, sizeof(int));

    /* Read command-line options. */
    while ((c = getopt(argc, argv, "s:i:rp:n:m:g")) != -1){
        switch (c){
        case 's':
            seed = atoi(optarg);
            break;
        case 'i':
            image_path = optarg;
            break;
        case 'r':
            gen_image = 1;
            break;
        case 'p':
            num_threads = atoi(optarg);
            break;
        case 'n':
            num_rows = strtol(optarg, 0, 10);
            break;
        case 'm':
            num_cols = strtol(optarg, 0, 10);
            break;
        case 'g':
            debug = 1;
            break;
        case '?':
            fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
            return -1;
        default:
            return -1;
        }
    }

    int* image = (int*)malloc(sizeof(int) * num_cols * num_rows);

    /* Seed such that we can always reproduce the same random vector */
    if (gen_image){
        srand(seed);
        generate_image(num_rows, num_cols, image);
    }
    else {
        read_image(image_path, num_rows, num_cols, image);
    }

    clock_gettime(CLOCK_MONOTONIC, &before);
    /* Do your thing here */

    histogram(histo, image, num_threads, num_rows * num_cols);

    /* Do your thing here */

    if (debug)
        print_histo(histo);

    clock_gettime(CLOCK_MONOTONIC, &after);
    double time = (double)(after.tv_sec - before.tv_sec) +
        (double)(after.tv_nsec - before.tv_nsec) / 1e9;

    printf("Histo took: % .6e seconds \n", time);
}