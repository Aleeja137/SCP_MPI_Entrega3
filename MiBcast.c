/*  SCP, Sistemas de Cómputo Paralelo -- GII (Ingeniería de Computadores)
    Entrega opcional - 1er ejercicio

    MiBcast.c

    Implementar una función para broadcast que crece logaritmicamente, 
    según un envío que se extiende por las coordenadas de un toro 2D

    Algoritmo:
    if (nodo_actual == origen)
        enviar_todos
    else
        switch (dirección_recv)
        case arriba:
            enviar_izq
            enviar_drc
            enviar_abajo (depende)
        case abajo:
            enviar_izq
            enviar_drc
            enviar_arriba (depende)
        case izq:
            enviar_drc (depende)
        case drc:
            enviar_izq (depende)
    
    Los que pone (depende) calcula según la distancia en la dimensión de envío respecto al nodo origen es menor al diámetro menos uno del toro.
    Cuando la cantidad de procesadores por dimensión es impar, o cuando se recibe el mensaje desde arriba o desde la izquierda,
    la distancia puede ser menor o igual al diámetro menos uno.
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define N 4000
#define REP 1

void MI_MPI_Bcast (int *mensaje, int ndat, int BCsrc, MPI_Comm C_toro)
{ 
    int pid_izq, pid_drc, pid_arriba, pid_abajo, coords[2], coords_src[2], pid;
    int  npr, dist = -1, diam = -1, parProcDim;
    int *mensajeRecv;

    // Cada proceso hace máximo 4 envíos
    MPI_Request requests[4];
    MPI_Status  stats[4];

    // Se calcula el diámetro, y si tenemos un número impar de procesadores por dimensión
    MPI_Comm_size (C_toro, &npr);
    diam = sqrt(npr)/2;
    if (diam*2 < sqrt(npr)) parProcDim = 1;
    else parProcDim = 0;

    // Se buscan los vecinos inmediatos en las cuatro direcciones
    MPI_Cart_shift(C_toro,1,1,&pid_izq,&pid_drc);
    MPI_Cart_shift(C_toro,0,1,&pid_arriba,&pid_abajo);

    // Se toman las coordenadas del nodo actual y la del nodo origen del mensaje
    MPI_Comm_rank(C_toro,&pid);
    MPI_Cart_coords(C_toro,pid,2,coords);  
    MPI_Cart_coords(C_toro,BCsrc,2,coords_src);

    // El nodo origen envía el mensaje en las cuatro direcciones siempre
    if (pid == BCsrc)
    {
        MPI_Isend(mensaje,ndat,MPI_INT,pid_arriba,0,C_toro,&requests[0]);
        MPI_Isend(mensaje,ndat,MPI_INT,pid_abajo,0, C_toro,&requests[1]);
        MPI_Isend(mensaje,ndat,MPI_INT,pid_izq,0,   C_toro,&requests[2]);
        MPI_Isend(mensaje,ndat,MPI_INT,pid_drc,0,   C_toro,&requests[3]);    
    }
    else
    {
        // Los demás nodos esperan el mensaje
        mensajeRecv = (int*) malloc (ndat*sizeof(int));
        MPI_Irecv(mensajeRecv,ndat,MPI_INT,MPI_ANY_SOURCE,0,C_toro,&requests[0]);
        MPI_Wait(&requests[0],&stats[0]);

        // Si el mensaje llega desde arriba, lo reenviamos a los lados 
        // Si estamos a menor distancia desde el origen que el diámetro, se manda abajo también
        if (stats[0].MPI_SOURCE == pid_arriba)
        {
            requests[0] = MPI_REQUEST_NULL; 
            dist = abs(coords[0] - coords_src[0]);
            if (dist > diam) dist = sqrt(npr) - dist;
            if (dist <= diam -1)  MPI_Isend(mensaje,ndat,MPI_INT,pid_abajo,0,C_toro,&requests[1]);
            else requests[1] = MPI_REQUEST_NULL; 
            MPI_Isend(mensaje,ndat,MPI_INT,pid_izq,0,   C_toro,&requests[2]);
            MPI_Isend(mensaje,ndat,MPI_INT,pid_drc,0,   C_toro,&requests[3]);
        }

        /* Si el mensaje llega desde abajo, lo reenviamos a los lados 
           Si estamos a menor distancia desde el origen que el diámetro 
           (la condición depende si hay un número par de procesadores por dimensión), 
           se manda arriba también */
        else if (stats[0].MPI_SOURCE == pid_abajo)
        {
            requests[1] = MPI_REQUEST_NULL;
            dist = abs(coords[0] - coords_src[0]);
            if (dist > diam) dist = sqrt(npr) - dist;
            if (dist < diam -1 || (parProcDim && (dist <= diam -1))) MPI_Isend(mensaje,ndat,MPI_INT,pid_arriba,0,C_toro,&requests[0]);
            else requests[0] = MPI_REQUEST_NULL;
            MPI_Isend(mensaje,ndat,MPI_INT,pid_izq,0,   C_toro,&requests[2]);
            MPI_Isend(mensaje,ndat,MPI_INT,pid_drc,0,   C_toro,&requests[3]);
        }

        // Si el mensaje llega desde la izquierda y 
        // estamos a menor distancia desde el origen que el diámetro, se manda a la derecha
        else if (stats[0].MPI_SOURCE == pid_izq)
        {
            requests[0] = requests[1] = requests[2] = MPI_REQUEST_NULL;
            dist = abs(coords[1] - coords_src[1]);
            if (dist > diam) dist = sqrt(npr) - dist;
            if (dist <= diam -1) MPI_Isend(mensaje,ndat,MPI_INT,pid_drc,0,C_toro,&requests[3]);
            else requests[3] = MPI_REQUEST_NULL; 
        }

        /* Si el mensaje llega desde la derecha, si estamos a menor distancia desde el origen que el diámetro 
           (la condición depende si hay un número par de procesadores por dimensión), 
           se manda a la izquierda también */
        else if (stats[0].MPI_SOURCE == pid_drc)
        {
            requests[0] = requests[1] = requests[3] = MPI_REQUEST_NULL;
            dist = abs(coords[1] - coords_src[1]);
            if (dist > diam) dist = sqrt(npr) - dist;
            if (dist < diam -1 || (parProcDim && (dist <= diam -1)))  MPI_Isend(mensaje,ndat,MPI_INT,pid_izq,0,C_toro,&requests[2]);
            else requests[2] = MPI_REQUEST_NULL;
        }        
    }

    // Se espera a que se hayan efectuado los mensajes
    MPI_Waitall(4,requests,stats);
    
} /*  MI_MPI_Bcast  */

void MI_MPI_Bcast_lineal (int *mensaje, int ndat, int BCsrc, MPI_Comm C_toro)
{ 
    int  npr, i, pid;
    int *mensajeRecv;

    // Necesitamos npr-1 envíos, pero por sencillez usamos npr requests
    MPI_Request *requests;
    MPI_Status  *stats;

    MPI_Comm_rank(C_toro,&pid);
    MPI_Comm_size (C_toro, &npr);
    
    if (pid == BCsrc)
    {
        requests = (MPI_Request*) malloc (npr * sizeof(MPI_Request));
        stats =    (MPI_Status*)  malloc (npr * sizeof(MPI_Status));

        for (i=0;i<npr;i++)
        {
            if (i != BCsrc) MPI_Isend(mensaje,ndat,MPI_INT,i,0,C_toro,&requests[i]);
            else requests[i] = MPI_REQUEST_NULL;
        }
        MPI_Waitall(npr,requests,stats);
    }
    else
    {
        requests = (MPI_Request*) malloc (sizeof(MPI_Request));
        stats =    (MPI_Status*)  malloc (sizeof(MPI_Status));
        mensajeRecv = (int*) malloc (ndat*sizeof(int));

        MPI_Irecv(mensajeRecv,ndat,MPI_INT,MPI_ANY_SOURCE,0,C_toro,requests);
        MPI_Wait(&requests[0],stats);
    }

} /*  MI_MPI_Bcast_lineal  */



int main (int argc, char *argv[]) 
{
    int pid, npr, BCsrc, i, A[N]; 
    int kdim;
    int ndims, dims[2], periods[2], reorder, pid_toro;
    double t0, t1, tex_mi_bcast = -1, tex_mi_bcast_lineal = -1, tex_bcast = -1;

    MPI_Comm C_toro;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &pid);
    MPI_Comm_size (MPI_COMM_WORLD, &npr);

    // Comprobar que npr es cuadrado
    kdim = sqrt(npr);
    if (kdim*kdim < npr || npr < 16)
    {
        if (pid == 0) printf ("\n OJO: número de procesos válidos: 16, 25, 36 ...\n\n\n");
        MPI_Finalize ();
        return (0);
    }

    // Crear la topología en forma de toro 2D
    ndims = 2;
    dims[0] = dims[1] = kdim;
    periods[0] = periods[1] = 1;
    reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,&C_toro);

    // Establecer el nodo origen e inicializar el mensaje a enviar
    BCsrc = 7;
    if (pid == BCsrc) 
    {
        for (i=0; i<N; i++) A[i] = pid; // N = 4000
        printf("\n BCsrc: %d \n",BCsrc);
    }

    // MI_BCast
    MPI_Barrier (MPI_COMM_WORLD);
    t0 = MPI_Wtime ();
    for (i=0; i<REP; i++) MI_MPI_Bcast (A, N, BCsrc, C_toro);
    MPI_Barrier (MPI_COMM_WORLD);
    t1 = MPI_Wtime ();
    tex_mi_bcast = (t1-t0) / REP;

    // MI_BCast_lineal
    MPI_Barrier (MPI_COMM_WORLD);
    t0 = MPI_Wtime ();
    for (i=0; i<REP; i++) MI_MPI_Bcast_lineal (A, N, BCsrc, C_toro);
    MPI_Barrier (MPI_COMM_WORLD);
    t1 = MPI_Wtime ();
    tex_mi_bcast_lineal = (t1-t0) / REP;

    // MPI_Bcast
    MPI_Barrier (MPI_COMM_WORLD);
    t0 = MPI_Wtime ();
    for (i=0; i<REP; i++) MPI_Bcast(A,N,MPI_INT,BCsrc,C_toro);
    MPI_Barrier (MPI_COMM_WORLD);
    t1 = MPI_Wtime ();
    tex_bcast = (t1-t0) / REP;

    if (pid == BCsrc) printf("\n MI_Bcast        tex: %f\n MI_Bcast_lineal tex: %f\n MPI_Bcast       tex: %f\n\n",tex_mi_bcast,tex_mi_bcast_lineal,tex_bcast);
    
    
    MPI_Finalize ();
    return (0);
} /*  main  */

