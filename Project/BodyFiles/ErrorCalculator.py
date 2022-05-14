import csv
from math import trunc

# CFG
reference_filepath = "out/out0.tsv"
approx_filepath = "out/bh_naive_out.tsv"

# LOGIC
def parse_bodies(filepath):
    bodies = []
    with open(filepath) as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        next(rd)
        for row in rd:
            body = [float(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5])]
            bodies.append(body)
    return bodies


ref_bodies = parse_bodies(reference_filepath)
approx_bodies = parse_bodies(approx_filepath)
ref_b_count = len(ref_bodies)
apx_b_count = len(approx_bodies)

if ref_b_count != apx_b_count:
    print(f"rip -> {ref_b_count} ref bodies VS {apx_b_count} approx bodies.")
    exit(1)

invalid_bodies = 0
max_p_error_body = max_v_error_body = -1
max_p_error_prcnt = max_v_error_prcnt = 0
avg_p_error_prcnt = avg_v_error_prcnt = 0

for i in range(ref_b_count):
    ref_body = ref_bodies[i]
    approx_body = approx_bodies[i]
    assert(ref_body[0] == approx_body[0])  # if masses are not the same: problem
    
    x_diff = ref_body[1] - approx_body[1]
    y_diff = ref_body[2] - approx_body[2]
    vel_x_diff = ref_body[3] - approx_body[3]
    vel_y_diff = ref_body[4] - approx_body[4]

    x_error_prcnt = abs(x_diff / ref_body[1]) * 100
    y_error_prcnt = abs(y_diff / ref_body[2]) * 100
    vel_x_error_prcnt = abs(vel_x_diff /  ref_body[3]) * 100
    vel_y_error_prcnt = abs(vel_y_diff /  ref_body[4]) * 100

    if x_error_prcnt > max_p_error_prcnt:
        max_p_error_prcnt = x_error_prcnt
        max_p_error_body = i

    if y_error_prcnt > max_p_error_prcnt:
        max_p_error_prcnt = y_error_prcnt
        max_p_error_body = i

    if vel_x_error_prcnt > max_v_error_prcnt:
        max_v_error_prcnt = vel_x_error_prcnt
        max_v_error_body = i

    if vel_y_error_prcnt > max_v_error_prcnt:
        max_v_error_prcnt = vel_y_error_prcnt
        max_v_error_body = i

    if x_error_prcnt > 1 or y_error_prcnt > 1 or vel_x_error_prcnt > 1 or vel_y_error_prcnt > 1:
        invalid_bodies += 1
    
    avg_p_error_prcnt += x_error_prcnt + y_error_prcnt
    avg_v_error_prcnt += vel_x_error_prcnt + vel_y_error_prcnt

avg_p_error_prcnt /= ref_b_count * 2
avg_v_error_prcnt /= ref_b_count * 2

avg_p_error_prcnt = round(avg_p_error_prcnt, 2)
avg_v_error_prcnt = round(avg_v_error_prcnt, 2)
max_p_error_prcnt = round(max_p_error_prcnt, 2)
max_v_error_prcnt = round(max_v_error_prcnt, 2)

print(f"\nInvalid bodies (error > 1% ): {invalid_bodies}\n")
print(f"AVG ERROR (%):\n  Position: {avg_p_error_prcnt}%\n  Velocity: {avg_v_error_prcnt}%\n")
print(f"MAX ERROR (%):\n  Position: {max_p_error_prcnt}% (body: {max_p_error_body})\n  Velocity: {max_v_error_prcnt}% (body: {max_v_error_body})\n")