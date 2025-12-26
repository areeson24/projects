***********************************************************************************************************
/*************** MACRO TO CREATE GROUPED MEANS +/- SEs LINE PLOT FOR REPEATED MEASURES DATA ***************/
***********************************************************************************************************;
/* Input dataset includes columns: &group | &x | &y | &se */
%macro plot_mean_se(
    data=,           /* Input dataset name */
    x=,              /* X-axis variable name (times) */
    y=,              /* Y-axis variable name (means) */
    se=,             /* Standard error variable to create error bars (leave blank to omit error bars) */
    group=,          /* Grouping variable (e.g. treatment) */
    group_order=,    /* Ordered list of groups for plotting separated by |: group1|group2|group3|etc. (optional) */
    colors=,         /* Colors for groups: color1 color2 color3 etc. (optional) */
    title=,          /* Plot title (optional) */
    subtitle=,       /* Plot subtitle (optional) */
    legend=,         /* Legend title (optional) */
    ylabel=,         /* Y-axis label (optional) */
    xlabel=,         /* X-axis label (optional) */
    ylim=,           /* Lower and upper plot bounds for y-axis: ymin ymax (optional) */
    ystep=,			 /* Step size for y-axis ticks (optional) */	
    outfile=         /* Output PDF file path + name (don't include .pdf extension, leave blank for no ODS output) */
);

/* Check variable names */
%let notfound=;
%let vars = &x &y &se &group;
%if %index(&data, .) > 0 %then %do;
	%let lib = %upcase(%scan(&data, 1, .));
	%let mem = %upcase(%scan(&data, 2, .));
%end;
%else %do;
	%let lib = WORK;
	%let mem = %upcase(&data);
%end;
%do i = 1 %to %sysfunc(countw(&vars));
    %let var = %scan(&vars, &i);
    proc sql noprint;
    	select count(*) into :exist_var trimmed
        from dictionary.columns
        where libname = "&lib"
          and memname = "&mem"
          and upcase(name) = upcase("&var");
    quit;
	%if &exist_var = 0 %then %let notfound = &notfound &var;
%end;
%if %length(&notfound) > 0 %then %do;
    %put ERROR: Variables not found in &data: &notfound.;
    %return;
%end;

/* Check group ordering */
%if %length(%sysfunc(strip(%superq(group_order)))) > 0 %then %do;
    proc sql noprint;
        select distinct strip(&group) into :data_groups separated by "|"
        from &data;
    quit;
    %let clean_group_order=;
    %let n_group_order = %sysfunc(countw(&group_order, |));
    %do i = 1 %to &n_group_order;
        %let this_group = %sysfunc(strip(%scan(&group_order, &i, |)));
        %let clean_group_order = &clean_group_order|&this_group;
    %end;
    %let clean_group_order = %substr(&clean_group_order, 2);
    %let n_data_groups = %sysfunc(countw(&data_groups, |));
    %let missing_groups=;
    %do i = 1 %to &n_data_groups;
        %let current = %sysfunc(strip(%scan(&data_groups, &i, |)));
        %if %index(|&clean_group_order.|, |&current.|) = 0 %then %do;
            %let missing_groups = &missing_groups &current;
        %end;
    %end;
    %if %length(&missing_groups) > 0 %then %do;
        %put ERROR: The following groups exist in the data but are missing from group_order:;
        %put &missing_groups;
        %return;
    %end;
%end;

/* Check colors */
%let n_colors = %sysfunc(countw(&colors));
%if &n_colors ^= &n_group_order %then %do;
	%put ERROR: The number of plotting colors (&n_colors) does not equal the number of groups (&n_group_order).;
	%return;
%end;

/* Check and extract y-axis bounds */
%let ymin=;
%let ymax=;
%if %length(%sysfunc(strip(&ylim))) > 0 %then %do;
    %let ymin = %scan(&ylim, 1, %str( ));
    %let ymax = %scan(&ylim, 2, %str( ));
    %if (%length(&ymax) = 0 or &ymax <= &ymin) %then %do;
    	%put ERROR: Y-axis limits are specified with 1 value missing or ymax <= ymin.;
    	%return;
    %end;
%end;

/* Subset and format data */
data tmp_df;
    set &data;
    length _group $100;
    _group = strip(&group);
    keep &x &y &se _group;
run;

/* Determine group order */
%if %length(%sysfunc(strip(%superq(group_order)))) > 0 %then %do;
    %let group_list=;
    %do i = 1 %to %sysfunc(countw(&group_order, |));
        %let group_clean = %sysfunc(strip(%scan(&group_order, &i, |)));
        %let group_list = &group_list|&group_clean;
    %end;
    %let group_list = %substr(&group_list, 2);
%end;
%else %do;
    proc sql noprint;
        select distinct _group into :group_list separated by "|"
        from tmp_df
        order by _group;
    quit;
%end;

%let n_groups = %sysfunc(countw(&group_list, |));

/* Plot error bounds? */
%if %length(%sysfunc(strip(%superq(se)))) > 0 %then %do;
    %let se_yn = Y;
    /* Create jitter map */
	data jitter_map;
	    length group_value $100 offset 8.;
	    %do i = 1 %to &n_groups;
	        group_value = "%scan(&group_list, &i, |)";
	        offset = %sysevalf((&i - (&n_groups + 1)/2) * 0.8);
	        output;
	    %end;
	run;
	/* Merge offsets and calculate error bounds */
	proc sql;
	    create table plot_df as
	    select a.*,
	           (a.&y - a.&se) as lower,
	           (a.&y + a.&se) as upper,
	           (a.&x + b.offset) as &x._jit
	    from tmp_df a
	    left join jitter_map b
	    on a._group = strip(b.group_value);
	quit;
%end;
%else %do;
	%let se_yn = N;
	data plot_df;
		set tmp_df;
	run;
%end;

/* Assign numeric group order based on group_list */
data plot_df;
    set plot_df;
    length _group_order 8;
    %do i = 1 %to &n_groups;
        if _group = "%scan(&group_list, &i, |)" then _group_order = &i;
    %end;
run;

/* Sort by group order and x */
proc sort data=plot_df;
    by _group_order &x;
run;

/* Get unique x values for axis */
proc sql noprint;
    select distinct &x into :x_list separated by " "
    from tmp_df
    order by &x;
quit;

/* ODS specs */
%let out = N;
%if %length(%sysfunc(strip(%superq(outfile)))) > 0 %then %do;
	%let out = Y;
	options nodate orientation=landscape;
	ods pdf file="&outfile..pdf" style=journal2 notoc;
%end;

/* Title and subtitle */
%if %length(%sysfunc(strip(%superq(title)))) > 0 %then %do;
    title "%superq(title)";
%end;
%else %do;
    title;
%end;
%if %length(%sysfunc(strip(%superq(subtitle)))) > 0 %then %do;
    title2 "%superq(subtitle)";
%end;
%else %do;
    title2;
%end;

/* Generate plot */
proc sgplot data=plot_df noautolegend;
	%if &se_yn = Y %then %do;
		scatter x=&x._jit y=&y / 
        	group=_group
        	yerrorlower=lower
        	yerrorupper=upper
        	markerattrs=(symbol=circlefilled size=3);
        series x=&x._jit y=&y / 
        	group=_group
        	markers
        	markerattrs=(symbol=circlefilled)
        	name="mean"
        	lineattrs=(pattern=solid thickness=1);
	%end;
	%else %if &se_yn = N %then %do;
		scatter x=&x y=&y / 
	        	group=_group
	        	markerattrs=(symbol=circlefilled size=3);
	    series x=&x y=&y / 
        	group=_group
        	markers
        	markerattrs=(symbol=circlefilled)
        	name="mean"
        	lineattrs=(pattern=solid thickness=1);
	%end;
    %if %length(&colors) > 0 %then %do;
        styleattrs datacontrastcolors=(&colors);
    %end;
    %if %length(%sysfunc(strip(%superq(xlabel)))) > 0 %then %do;
    	xaxis 
    		label="%superq(xlabel)"
    		offsetmin=0.07 offsetmax=0.07
          	values=(&x_list);
	%end;
	%else %do;
    	xaxis 
    		display=(nolabel)
            offsetmin=0.07 offsetmax=0.07
            values=(&x_list);
	%end;
	%if %length(%sysfunc(strip(%superq(ylabel)))) > 0 %then %do;
        yaxis 
            label="%superq(ylabel)"
            %if %length(&ylim) > 0 %then %do;
                %if %length(&ystep) > 0 %then %do;
                    min=&ymin max=&ymax values=(&ymin TO &ymax BY &ystep)
                %end;
                %else %do;
                    min=&ymin max=&ymax
                %end;
            %end;;
    %end;
    %else %do;
        yaxis 
            display=(nolabel)
            %if %length(&ylim) > 0 %then %do;
                %if %length(&ystep) > 0 %then %do;
                    min=&ymin max=&ymax values=(&ymin TO &ymax BY &ystep)
                %end;
                %else %do;
                    min=&ymin max=&ymax
                %end;
            %end;;
    %end;
    keylegend "mean" / title="%superq(legend)";
run;

title;
title2;

%if &out = Y %then %do;
	ods pdf close;
%end;

%mend;

********************************************************************************************;

/* TEST */

libname tmp "/home/u64253485/datasets/";

data df;
	set tmp.f14_2_1_1_1a;
run;

%plot_mean_se(
    data=df, 
    x=week, 
    y=mean,
    se=se,
    group=trt, 
    group_order=Placebo TID|Pirfenidone TID|LYT-100 550 mg TID|LYT-100 825 mg TID,
    colors=black gray steelblue green,
    title=Observed Mean FVC Changes Over Time - Part A, 
    subtitle=Full Analysis Set, 
    legend=Treatment, 
    xlabel=Week, 
    ylabel=Mean (SE) Change in FVC (mL) from Baseline,
    ylim=-120 60,
    ystep=20,
    outfile=/home/u64253485/test
);

