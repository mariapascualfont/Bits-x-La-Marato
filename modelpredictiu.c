#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_LINE 1024
#define MAX_KIR_TYPES 100
#define MAX_SUBTYPES 3000

typedef struct {
    char name[100];
    char main_type[50];
    int count;
    int relapse_cases;
    double relapse_prob;
    double frequency_percent;
} AlleleData;

AlleleData alleles[MAX_SUBTYPES];
AlleleData main_types[MAX_KIR_TYPES];
int total_alleles = 0;
int unique_subtypes = 0;
int unique_main_types = 0;

void extract_main_type(const char *allele_name, char *main_type) {
    int i = 0, j = 0;
    while (allele_name[i] != '*' && allele_name[i] != '\0') {
        main_type[j++] = allele_name[i++];
    }
    main_type[j] = '\0';
}

void process_file() {
    FILE *file = fopen("Allelelist.txt", "r");
    if (!file) {
        printf("Error: No se pudo abrir Allelelist.txt\n");
        exit(1);
    }

    char line[MAX_LINE];
    char allele_id[50], allele_name[100];

    // Saltar encabezado
    fgets(line, sizeof(line), file);

    while (fgets(line, sizeof(line), file)) {
        if (sscanf(line, "%[^,],%s", allele_id, allele_name) == 2) {
            total_alleles++;

            // Buscar subtipo existente
            int found = 0;
            for (int i = 0; i < unique_subtypes; i++) {
                if (strcmp(alleles[i].name, allele_name) == 0) {
                    alleles[i].count++;
                    found = 1;
                    break;
                }
            }

            // Si no existe, agregar nuevo
            if (!found && unique_subtypes < MAX_SUBTYPES) {
                strcpy(alleles[unique_subtypes].name, allele_name);
                alleles[unique_subtypes].count = 1;
                alleles[unique_subtypes].relapse_cases = rand() % 6; // Simulado
                extract_main_type(allele_name, alleles[unique_subtypes].main_type);
                unique_subtypes++;
            }

            // Actualizar tipo principal
            char main_type[50];
            extract_main_type(allele_name, main_type);
            found = 0;
            for (int i = 0; i < unique_main_types; i++) {
                if (strcmp(main_types[i].main_type, main_type) == 0) {
                    main_types[i].count++;
                    // Sumar recaídas del subtipo correspondiente
                    for (int j = 0; j < unique_subtypes; j++) {
                        if (strcmp(alleles[j].name, allele_name) == 0) {
                            main_types[i].relapse_cases += alleles[j].relapse_cases;
                            break;
                        }
                    }
                    found = 1;
                    break;
                }
            }

            if (!found && unique_main_types < MAX_KIR_TYPES) {
                strcpy(main_types[unique_main_types].main_type, main_type);
                strcpy(main_types[unique_main_types].name, main_type);
                main_types[unique_main_types].count = 1;
                for (int j = 0; j < unique_subtypes; j++) {
                    if (strcmp(alleles[j].name, allele_name) == 0) {
                        main_types[unique_main_types].relapse_cases = alleles[j].relapse_cases;
                        break;
                    }
                }
                unique_main_types++;
            }
        }
    }
    fclose(file);
}

void calculate_stats() {
    for (int i = 0; i < unique_subtypes; i++) {
        alleles[i].frequency_percent = (double)alleles[i].count / total_alleles * 100.0;
        alleles[i].relapse_prob = (double)alleles[i].relapse_cases / alleles[i].count * 100.0;
    }

    for (int i = 0; i < unique_main_types; i++) {
        main_types[i].frequency_percent = (double)main_types[i].count / total_alleles * 100.0;
        if (main_types[i].count > 0) {
            main_types[i].relapse_prob = (double)main_types[i].relapse_cases / main_types[i].count * 100.0;
        }
    }
}

void sort_by_probability(AlleleData *data, int count) {
    for (int i = 0; i < count - 1; i++) {
        for (int j = i + 1; j < count; j++) {
            if (data[j].relapse_prob > data[i].relapse_prob) {
                AlleleData temp = data[i];
                data[i] = data[j];
                data[j] = temp;
            }
        }
    }
}

void generate_csv(const char *filename, AlleleData *data, int count, int is_main_type) {
    FILE *csv = fopen(filename, "w");
    if (!csv) {
        printf("Error creando %s\n", filename);
        return;
    }

    if (is_main_type) {
        fprintf(csv, "Main_Type,Total_Count,Relapse_Cases,Probability(%%),Frequency(%%)\n");
        for (int i = 0; i < count; i++) {
            fprintf(csv, "%s,%d,%d,%.2f,%.4f\n",
                   data[i].main_type,
                   data[i].count,
                   data[i].relapse_cases,
                   data[i].relapse_prob,
                   data[i].frequency_percent);
        }
    } else {
        fprintf(csv, "Allele,Main_Type,Count,Relapse_Cases,Probability(%%),Frequency(%%)\n");
        for (int i = 0; i < count; i++) {
            fprintf(csv, "%s,%s,%d,%d,%.2f,%.4f\n",
                   data[i].name,
                   data[i].main_type,
                   data[i].count,
                   data[i].relapse_cases,
                   data[i].relapse_prob,
                   data[i].frequency_percent);
        }
    }
    fclose(csv);
    printf("Archivo generado: %s (%d registros)\n", filename, count);
}

void display_all_results() {
    printf("\n=======================================================\n");
    printf("TODOS LOS RESULTADOS - %d SUBTIPOS ANALIZADOS\n", unique_subtypes);
    printf("=======================================================\n");

    printf("\n1. TODOS LOS SUBTIPOS ordenados por probabilidad (mayor a menor):\n");
    printf("=================================================================\n");
    printf("%-35s %-12s %-8s %-12s %-10s\n", 
           "Allelo", "Tipo KIR", "Count", "Prob(%)", "Freq(%)");
    printf("=================================================================\n");

    for (int i = 0; i < unique_subtypes; i++) {
        printf("%-35s %-12s %-8d %-12.2f %-10.4f\n",
               alleles[i].name,
               alleles[i].main_type,
               alleles[i].count,
               alleles[i].relapse_prob,
               alleles[i].frequency_percent);
    }

    printf("\n\n2. TODOS LOS TIPOS PRINCIPALES ordenados por probabilidad:\n");
    printf("=============================================================\n");
    printf("%-15s %-10s %-12s %-10s %-10s\n", 
           "Tipo KIR", "Total", "Recaídas", "Prob(%)", "Freq(%)");
    printf("=============================================================\n");

    for (int i = 0; i < unique_main_types; i++) {
        printf("%-15s %-10d %-12d %-10.2f %-10.4f\n",
               main_types[i].main_type,
               main_types[i].count,
               main_types[i].relapse_cases,
               main_types[i].relapse_prob,
               main_types[i].frequency_percent);
    }

    printf("\n\nRESUMEN ESTADÍSTICO:\n");
    printf("===================\n");
    printf("Total alelos procesados: %d\n", total_alleles);
    printf("Subtipo únicos: %d\n", unique_subtypes);
    printf("Tipos KIR principales: %d\n", unique_main_types);

    // Estadísticas adicionales
    double avg_prob = 0, max_prob = 0, min_prob = 100;
    char max_kir[50], min_kir[50];

    for (int i = 0; i < unique_main_types; i++) {
        avg_prob += main_types[i].relapse_prob;
        if (main_types[i].relapse_prob > max_prob) {
            max_prob = main_types[i].relapse_prob;
            strcpy(max_kir, main_types[i].main_type);
        }
        if (main_types[i].relapse_prob < min_prob) {
            min_prob = main_types[i].relapse_prob;
            strcpy(min_kir, main_types[i].main_type);
        }
    }
    avg_prob /= unique_main_types;

    printf("Probabilidad promedio: %.2f%%\n", avg_prob);
    printf("Mayor probabilidad: %s (%.2f%%)\n", max_kir, max_prob);
    printf("Menor probabilidad: %s (%.2f%%)\n", min_kir, min_prob);
}

int main() {
    printf("Procesando Allelelist.txt...\n");
    process_file();

    printf("Calculando estadísticas...\n");
    calculate_stats();

    printf("Ordenando por probabilidad...\n");
    sort_by_probability(alleles, unique_subtypes);
    sort_by_probability(main_types, unique_main_types);

    printf("\nGenerando archivos CSV...\n");
    generate_csv("kir_ALL_subtypes_by_probability.csv", alleles, unique_subtypes, 0);
    generate_csv("kir_ALL_main_types_by_probability.csv", main_types, unique_main_types, 1);

    printf("\n¿Desea ver TODOS los resultados en pantalla? (s/n): ");
    char respuesta;
    scanf(" %c", &respuesta);

    if (respuesta == 's' || respuesta == 'S') {
        display_all_results();
    } else {
        printf("\nPuede abrir los archivos CSV en Excel para ver todos los resultados.\n");
    }

    printf("\nArchivos generados:\n");
    printf("1. kir_ALL_subtypes_by_probability.csv - Todos los subtipos (%d)\n", unique_subtypes);
    printf("2. kir_ALL_main_types_by_probability.csv - Todos los tipos principales (%d)\n", unique_main_types);
    printf("\nAbra estos archivos en Excel para análisis completo.\n");

    return 0;
}