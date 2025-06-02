import numpy as np
import matplotlib.pyplot as plt

# --- 精确解函数 ---
def exact_solution(x, t):
    return np.sin(np.pi * x) * np.exp(-np.pi**2 * t) + x

# --- 近似解函数 ---
def approximate_solution(x, t, coeffs_a, N):
    """
    计算近似解 theta_approx = sin(pi*x) + x + sum_{j=1}^{N} a_j(t) (x^j - x^{j+1})
    coeffs_a: 包含 a_1(t), a_2(t), ..., a_N(t) 的列表或数组
    N: 系数的个数
    """
    if N == 0: # 对应原始数据中 N=1 的情况，coeffs_a 只有一个元素
        if len(coeffs_a) == 0: # 处理空coeffs_a的情况
             return np.sin(np.pi * x) + x
        term_sum = coeffs_a[0] * (x**1 - x**2)
    else:
        term_sum = 0.0
        for j in range(N): # j 从 0 到 N-1
            # coeffs_a[j] 对应 a_{j+1}(t)
            # (x^{j+1} - x^{j+2})
            term_sum += coeffs_a[j] * (x**(j + 1) - x**(j + 2))
    return np.sin(np.pi * x) + x + term_sum


# --- 均方根误差计算函数 ---
def calculate_rmse(numerical_sol, exact_sol):
    return np.sqrt(np.mean((numerical_sol - exact_sol)**2))

# --- 主程序 ---
def main():
    filepath = '../output.txt' # 数据文件路径
    time_to_evaluate = 0.20 # 所有数据都对应这个时间点

    # 读取数据
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"错误: 文件 '{filepath}' 未找到。请确保文件路径正确。")
        return

    data_by_method = {
        "Galerkin": [],
        "Subdomain": [],
        "Collocation": []
    }

    # 解析数据
    # 假设每种方法有4行数据，对应 N=1, 3, 5, 7
    # 并且每行第一个是时间，后面是 a_j 系数
    n_values_actual = [1, 3, 5, 7] # 实际的N值，用于近似解函数
    # 注意：原始数据中N=1意味着 coeffs_a 只有一个元素，对应 a_1
    # N=3 意味着 coeffs_a 有三个元素，对应 a_1, a_2, a_3，以此类推

    method_order = ["Galerkin", "Subdomain", "Collocation"]
    current_method_idx = 0
    line_idx_in_method = 0

    parsed_data = []
    for line_str in lines:
        if not line_str.strip(): # 跳过空行
            continue
        parts = list(map(float, line_str.strip().split()))
        time_val = parts[0]
        coeffs = parts[1:]
        parsed_data.append({'time': time_val, 'coeffs': coeffs, 'N_actual': n_values_actual[line_idx_in_method]})

        line_idx_in_method += 1
        if line_idx_in_method == len(n_values_actual):
            data_by_method[method_order[current_method_idx]] = parsed_data
            parsed_data = []
            line_idx_in_method = 0
            current_method_idx += 1
            if current_method_idx >= len(method_order):
                break # 所有方法数据读取完毕

    if current_method_idx < len(method_order) and parsed_data: # 处理最后一种方法的数据（如果文件行数不足）
         data_by_method[method_order[current_method_idx]] = parsed_data


    # --- 绘图 ---
    fig, axes = plt.subplots(nrows=len(n_values_actual), ncols=len(method_order),
                             figsize=(12, 18), sharex=True, sharey=True)
    fig.suptitle(f'Numerical vs. Exact Solution at t = {time_to_evaluate:.2f}', fontsize=16)

    x_points = np.linspace(0, 1, 100) # 用于绘图和误差计算的空间点

    for col_idx, method_name in enumerate(method_order):
        method_data_list = data_by_method.get(method_name, [])
        if not method_data_list:
            print(f"警告: 没有找到方法 '{method_name}' 的数据。")
            for row_idx in range(len(n_values_actual)):
                ax = axes[row_idx, col_idx] if len(n_values_actual) > 1 else axes[col_idx]
                ax.text(0.5, 0.5, "No Data", ha='center', va='center', fontsize=10, color='red')
                if row_idx == len(n_values_actual) -1 :
                     ax.set_xlabel('x')
                if col_idx == 0:
                    ax.set_ylabel(f'N={n_values_actual[row_idx]}\nθ(x,t)')
                if row_idx == 0:
                    ax.set_title(method_name)
            continue


        for row_idx, N_actual in enumerate(n_values_actual):
            ax = axes[row_idx, col_idx] if len(n_values_actual) > 1 else axes[col_idx]

            # 查找对应N的数据
            data_entry = None
            for entry in method_data_list:
                # 我们的数据是按顺序读取的，第row_idx个数据对应当前的N_actual
                if len(method_data_list) > row_idx : #确保列表有足够的元素
                    data_entry = method_data_list[row_idx]
                    break
            
            if data_entry is None or data_entry['N_actual'] != N_actual:
                print(f"警告: 方法 '{method_name}' 中 N={N_actual} 的数据缺失或不匹配。")
                ax.text(0.5, 0.5, f"Data Missing\nfor N={N_actual}", ha='center', va='center', fontsize=9, color='red')
                if row_idx == len(n_values_actual) -1 :
                     ax.set_xlabel('x')
                if col_idx == 0:
                    ax.set_ylabel(f'N={N_actual}\nθ(x,t)')
                if row_idx == 0:
                    ax.set_title(method_name)
                continue


            coeffs_a = data_entry['coeffs']
            
            # 计算数值解和精确解
            numerical_sol_values = np.array([approximate_solution(x, time_to_evaluate, coeffs_a, N_actual) for x in x_points])
            exact_sol_values = exact_solution(x_points, time_to_evaluate)

            # 绘图
            ax.plot(x_points, exact_sol_values, 'k-', label='Exact')
            ax.plot(x_points, numerical_sol_values, 'r--', label='Numerical')
            
            # 计算并显示RMSE
            rmse = calculate_rmse(numerical_sol_values, exact_sol_values)
            ax.text(0.05, 0.05, f'RMSE: {rmse:.2e}', transform=ax.transAxes, fontsize=8,
                    bbox=dict(boxstyle='round,pad=0.3', fc='wheat', alpha=0.5))

            if row_idx == 0:
                ax.set_title(method_name)
            if col_idx == 0:
                ax.set_ylabel(f'N={N_actual}\nθ(x,t)',fontsize = 14)
            if row_idx == len(n_values_actual) -1 :
                ax.set_xlabel('x')

            ax.grid(True, linestyle=':', alpha=0.7)
            if row_idx == 0 and col_idx == len(method_order) -1: # 只在一个子图上显示图例
                ax.legend(fontsize=8, loc='upper right')


    plt.tight_layout(rect=[0, 0, 1, 0.96]) # 调整布局以适应主标题
    plt.savefig('solution_comparison_rmse.png')
    plt.show()

if __name__ == '__main__':
    main()