function output_config(filepath, dataset, remove_thd, persist_thd, keep_top)
fp = fopen([filepath int2str(dataset) '.txt'], 'w');
fprintf(fp, 'Dataset #: %d\n', dataset);
fprintf(fp, 'Threshold for reducing input %f\n', remove_thd);
fprintf(fp, 'Threshold for persistence simplification %f\n', persist_thd);
fprintf(fp, 'Keep to %d persistence branches.\n', keep_top);