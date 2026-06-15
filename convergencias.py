import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv

# ── Leer CSV ────────────────────────────────────────────────────────────────
def leer_csv(path):
    iters, makespans = [], []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            iters.append(int(row['Iteracion']))
            makespans.append(int(row['Makespan']))
    return iters, makespans

iters_p, ms_p = leer_csv('convergencia_Pequena_5x4.csv')
iters_m, ms_m = leer_csv('convergencia_Mediana_10x6.csv')
iters_g, ms_g = leer_csv('convergencia_Grande_20x5.csv')

# ── Configuración visual ─────────────────────────────────────────────────────
COLOR_LINEA  = '#1a56db'   # azul para la curva
COLOR_INICIO = '#e3452f'   # rojo para el punto inicial
COLOR_MEJOR  = '#16a34a'   # verde para el mejor valor encontrado
FONDO_PANEL  = '#f8fafc'   # gris muy claro para el fondo de cada subplot

# ── Figura con 3 subplots ────────────────────────────────────────────────────
# GridSpec: 2 filas x 2 columnas
# - Pequeña  → fila 0, col 0
# - Mediana  → fila 0, col 1
# - Grande   → fila 1, cols 0-1 (ocupa toda la fila, es la más interesante)
fig = plt.figure(figsize=(15, 12), facecolor='white')
fig.suptitle(
    'Curvas de Convergencia — Fireworks Algorithm (PFSP)\n'
    'Grupo 4 | Análisis y Diseño de Algoritmos',
    fontsize=14, fontweight='bold', y=0.97, color='#1e293b'
)

gs = gridspec.GridSpec(
    2, 2, figure=fig,
    hspace=0.42, wspace=0.32,
    left=0.07, right=0.97, top=0.91, bottom=0.07
)

# ── Datos por instancia ──────────────────────────────────────────────────────
instancias = [
    {
        'iters': iters_p,
        'ms':    ms_p,
        'titulo':  'Instancia Pequeña (5 trabajos × 4 máquinas)',
        'archivo': 'instancia1_bas1.txt',
        'pos':     gs[0, 0],   # subplot superior izquierdo
    },
    {
        'iters': iters_m,
        'ms':    ms_m,
        'titulo':  'Instancia Mediana (10 trabajos × 6 máquinas)',
        'archivo': 'instancia2_car5.txt',
        'pos':     gs[0, 1],   # subplot superior derecho
    },
    {
        'iters': iters_g,
        'ms':    ms_g,
        'titulo':  'Instancia Grande (20 trabajos × 5 máquinas)',
        'archivo': 'instancia3_reC01.txt',
        'pos':     gs[1, :],   # subplot inferior, ocupa las 2 columnas
    },
]

# ── Dibujar cada subplot ─────────────────────────────────────────────────────
for inst in instancias:
    iters = inst['iters']
    ms    = inst['ms']
    ax    = fig.add_subplot(inst['pos'])

    # Fondo y grilla
    ax.set_facecolor(FONDO_PANEL)
    ax.grid(True, linestyle='--', alpha=0.5, color='#cbd5e1', zorder=0)
    for spine in ax.spines.values():
        spine.set_edgecolor('#94a3b8')
        spine.set_linewidth(0.8)

    # Línea escalonada (step='post': el valor cambia al inicio del intervalo)
    ax.step(iters, ms, where='post',
            color=COLOR_LINEA, linewidth=1.8, zorder=3, label='Mejor makespan')

    # Punto rojo = valor inicial (primera iteración)
    ax.scatter(iters[0], ms[0],
               color=COLOR_INICIO, s=60, zorder=5,
               label=f'Inicial: {ms[0]}')

    # Punto verde estrella = mejor valor encontrado en toda la ejecución
    best_ms   = min(ms)
    best_iter = next(i for i, v in zip(iters, ms) if v == best_ms)
    ax.scatter(best_iter, best_ms,
               color=COLOR_MEJOR, s=70, marker='*', zorder=5,
               label=f'Mejor: {best_ms} (iter {best_iter})')

    # Línea horizontal punteada en el mejor valor (referencia visual)
    ax.axhline(y=best_ms, color=COLOR_MEJOR,
               linestyle=':', linewidth=1.0, alpha=0.7, zorder=2)

    # Títulos y etiquetas
    ax.set_title(inst['titulo'],
                 fontsize=10.5, fontweight='bold', color='#1e293b', pad=7)
    ax.set_xlabel(f'Iteración  |  Archivo: {inst["archivo"]}',
                  fontsize=8.5, color='#475569')
    ax.set_ylabel('Mejor Makespan', fontsize=8.5, color='#475569')
    ax.tick_params(labelsize=8, colors='#475569')

    # Leyenda
    leg = ax.legend(fontsize=8, loc='upper right',
                    framealpha=0.9, edgecolor='#cbd5e1', fancybox=False)
    leg.get_frame().set_linewidth(0.8)

    # Padding extra en Y para instancias con curva plana (evita que quede aplastada)
    ymin, ymax = ax.get_ylim()
    if (ymax - ymin) < 5:
        ax.set_ylim(ymin - 10, ymax + 10)

# ── Guardar ──────────────────────────────────────────────────────────────────
plt.savefig('curvas_convergencia.png', dpi=180, bbox_inches='tight', facecolor='white')
print("Imagen guardada en: curvas_convergencia.png")
