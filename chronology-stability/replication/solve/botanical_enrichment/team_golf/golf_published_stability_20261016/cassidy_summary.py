"""Compact diagnostic readout of derivative-based H continuation."""
from pathlib import Path
import hashlib
import json
import numpy as np
from threadpoolctl import threadpool_limits
import cassidy_h_continuation as study

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    result_path = HERE/'cassidy_h_results.json'
    result = json.loads(result_path.read_text())
    assert result['status']=='complete'
    fd = result['finite_difference_checks']
    out = dict(status='complete', source_sha256=sha(result_path), rendering_or_summary_source_sha256=sha(__file__),
        derivative_checks=dict(number=len(fd), maximum_relative_error=max(r['derivative_relative_error'] for r in fd),
            maximum_mixed_derivative_error=max(r['mixed_derivative_max_error'] for r in fd),
            maximum_hessian_action_error=max(r['hessian_action_max_error'] for r in fd)),
        event_checks=result['aggregate'], counts=result['counts'],
        half_nudge_checks=dict(number=len(result['path_refinement_checks']),
            maximum_node_rms=max(r['nudge_half_endpoint_rms'] for r in result['path_refinement_checks']),
            maximum_objective_difference=max(abs(r['objective_difference']) for r in result['path_refinement_checks']),
            order_disagreements=sum(not r['same_relation'] for r in result['path_refinement_checks'])),
        planes=[], baseline_derivative_sensitivity=[])
    for row in result['settings']:
        cells=row['cells']
        z=np.array([c['combined_herbal_preference']['gap_percent'] for c in cells])
        old=np.array([c['original_herbal_preference']['gap_percent'] for c in cells])
        disagreement=sum(len(set(np.sign(p['herbal_preference']['gap_percent']) for p in c['paths']))>1 for c in cells)
        disagreeing_indices=[i for i,c in enumerate(cells)
            if len(set(np.sign(p['herbal_preference']['gap_percent']) for p in c['paths']))>1]
        native=np.asarray(result['reference_coords3'])[result['h_index'],row['plane']]
        example=None
        if disagreeing_indices:
            near=min(disagreeing_indices,key=lambda i:
                (row['xs'][i%25]-native[0])**2+(row['ys'][i//25]-native[1])**2)
            c=cells[near]
            example=dict(cell_index=near,coordinates=[row['xs'][near%25],row['ys'][near//25]],
                original_gap_percent=c['original_herbal_preference']['gap_percent'],
                combined_gap_percent=c['combined_herbal_preference']['gap_percent'],
                path_gaps_percent={p['name']:p['herbal_preference']['gap_percent'] for p in c['paths']},
                node_rms_across_paths=c['maximum_path_node_rms'])
        out['planes'].append(dict(plane=row['plane'],cells=len(cells),positive3a_first=int(np.sum(z>1e-9)),
            negative3b_first=int(np.sum(z < -1e-9)), ties=int(np.sum(abs(z)<=1e-9)),
            path_preference_disagreement_cells=disagreement,
            original_map_preference_changes=int(np.sum(np.sign(old)!=np.sign(z))),
            maximum_gap_change_from_original_percent=float(np.max(abs(old-z))),
            path_any_seed_order_disagreement_cells=sum(c['path_order_changes']>0 for c in cells),
            path_any_seed_nodes_different_cells=sum(c['maximum_path_node_rms']>1e-8 for c in cells),
            return_any_seed_order_change_cells=sum(c['return_order_changes']>0 for c in cells),
            return_any_seed_nodes_different_cells=sum(c['maximum_return_node_rms']>1e-8 for c in cells),
            maximum_path_node_rms=max(c['maximum_path_node_rms'] for c in cells),
            maximum_return_node_rms=max(c['maximum_return_node_rms'] for c in cells),
            native_H_gap_percent=cells[row['baseline_cell_index']]['combined_herbal_preference']['gap_percent'],
            nearest_native_H_path_preference_disagreement=example,
            terminal4c_first_cells=sum(c['combined_terminal_preference']['gap_percent']<0 for c in cells)))
    # This small exact derivative calculation makes no new curve fits.
    data=study.grid.COMMON.load_data(); x,_=study.grid.COMMON.embed(data['pool_ds']['combined'])
    old=json.loads((study.OLD/'expanded_influence_grid.json').read_text())
    basins=json.loads((HERE.parent/'golf_stability_20261016/basins.json').read_text())
    base=next(s for s in basins['settings'] if s['cohort']=='pool50' and s['lambda']==.01)
    enriched=next(s for s in old['seed_sets'] if s['penalty']==.01)
    seeds=[(r['start'],np.asarray(r['nodes'])) for r in base['baseline_runs']]
    seeds += [(r['name'],np.asarray(r['nodes'])) for r in enriched['added_full_models']]
    with threadpool_limits(limits=1):
        for name,nodes in seeds:
            coeff,assigned,eigen=study.derivative(x,nodes,result['h_index'])
            out['baseline_derivative_sensitivity'].append(dict(seed=name,
                rms_node_displacement_per_unit_H_displacement=float(np.sqrt(np.mean(coeff*coeff))),
                maximum_node_displacement_per_unit_H_displacement=float(np.max(abs(coeff))),
                minimum_hessian_eigenvalue=eigen,H_assigned_node=int(assigned[result['h_index']]),
                note='Local spatial sensitivity within the current node-assignment region; not a time exponent.'))
    study.save('cassidy_h_summary.json',out)
    print(json.dumps({k:v for k,v in out.items() if k!='baseline_derivative_sensitivity'},indent=2))


if __name__=='__main__':
    main()
