#' # functions for working wiht eugene's sequence databases

#' Return the seuqnece entry corresponding to a given hi name (v$s[[idx]]$h)
match_hi_name = function(name, ...){

  dbs = list(...)

  db_idx=1
  while(db_idx <= length(dbs)){
    db = dbs[[ db_idx ]]; db_idx = db_idx+1;

    search_names = unlist(lapply(h3n2,
                                 function(v){
                                   lapply(v$s,
                                          function(seq_entry){
                                            seq_entry$h
                                          })
                                 }
    )
    )

    name_match_idx1 = which(unlist(lapply(search_names, function(ns)isTRUE(name %in% unlist(ns) ))))

    name_match_idx2 = lapply(name_match_idx1,
                             function(idx1){
                               which(unlist(search_names[[idx1]] == name))[[1]]
                             }
    )





    for (i in seq_along(name_match_idx1)){
      potential_match = db[[name_match_idx1[[i]]]]
      if (potential_match$N == name){

        if (!(name == potential_match$N & hash == potential_match$s[[ hash_match_idx2[[i]]  ]]$H)) stop('What\'s going on here...')

        return(potential_match$s[[ hash_match_idx2[[i]] ]])
      }
    }

  }

  warning('No match found')

  return(NULL)





}







#' Match many hashes
#'
#' iterates over database to find seq_entries corresponding to hashes
#'@export
match_many_hashes = function(Hs, db){

  out_seq_entries = lapply(Hs, function(h) 'missing')

  for (v in db){
    for (s in v$s){
      if (!is.null(s$H)){
        idxs = which(Hs == s$H)
        if (length(idxs > 0)){
          for (i in idxs) out_seq_entries[i] = list(s)
        }

      }
    }
  }

  n_missing = sum(sapply(out_seq_entries, function(e)!is.list(e)))
  if (n_missing > 0) warning('Could not find match for ', n_missing, ' hashes.')

  out_seq_entries
}



#' Find sequence entry with matching R (hash and name)
#'@export
match_ref = function(R, ..., ignore_name = F, return_full_virus = F){

  dbs = list(...)

  hash = R$H
  name = R$N

  db_idx=1
  while(db_idx <= length(dbs)){
    db = dbs[[ db_idx ]]; db_idx = db_idx+1;

    search_hashes = lapply(db, function(v)lapply(v$s, function(x)return(x$H)))

    hash_match_idx1 = which(unlist(lapply(search_hashes, function(hs)isTRUE(hash %in% unlist(hs) ))))

    hash_match_idx2 = lapply(hash_match_idx1,
                             function(idx1){
                               which(unlist(search_hashes[[idx1]] == hash))[[1]]
                             }
    )



    for (i in seq_along(hash_match_idx1)){
      potential_match = db[[hash_match_idx1[[i]]]]
      if (potential_match$N == name | ignore_name){

        if (!( (name == potential_match$N | ignore_name) & hash == potential_match$s[[ hash_match_idx2[[i]]  ]]$H)) stop('What\'s going on here...')

        if (return_full_virus) return(potential_match)
        else return(potential_match$s[[ hash_match_idx2[[i]] ]])
      }
    }

  }

  warning('No match found for: name = ', name, '  hash = ', hash)

  return(NULL)
}

#' Make sure seqquence entry has sequence
#'@export
with_sequence.seq_entry = function(seq_entry, ..., ignore_name = F){
  if (!is.null(seq_entry$a)) return(seq_entry)

  seq_match = match_ref(seq_entry$R, ..., ignore_name = ignore_name)

  if (is.null(seq_match)){
    warning('No match found for hash: ', seq_entry$R$H)
  }

  seq_entry[c('a', 's', 'n', 't', 'c')] = seq_match[c('a', 's', 'n', 't', 'c')]


  seq_entry
}

#' Make sure all sequence_entries for virus have sequence
#'@export
with_sequence = function(v, ..., ignore_name = F){

  for (s_idx in seq_along(v$s)){

    if (is.null(v$s[[s_idx]]$a)){


      v$s[[s_idx]] = with_sequence.seq_entry(v$s[[s_idx]], ..., ignore_name = ignore_name)

    }
  }
  v
}

#' Get trimmed sequence from sequence entry
#'@export
trim_sequence = function(seq_entry){

  if (is.null(seq_entry$s)){
    return(seq_entry$a)
  }
  return(stringr::str_sub(seq_entry$a, start = 1-(seq_entry$s)))
}

# return list with each trimmed sequence and clade for a virus
#'@export
get_sequences <- function(v){
  seqlist = list()

  for (s_idx in seq_along(v$s)){
    seq_entry = v$s[[s_idx]]
    if (!is.null(seq_entry$a)){
      trimmed_seq = trim_sequence(seq_entry)
      seqlist[[s_idx]] = list(sequence = trimmed_seq, clade = seq_entry$c)
    }
    else{
      stop('null sequence!')
    }

  }
  seqlist
}




# tree stuff
###########

#' Get list of all virus info in tree
#'@export
flatten_eutree = function(eutree){

  local_tree_vs = list()


  if ('t' %in% names(eutree)){

    for (entry in eutree$t ){

      if ('n' %in% names(entry)){

        if ('t' %in% names(entry)) stop('wtf')

        local_tree_vs = c(local_tree_vs, list(entry))
      }
      else{
        local_tree_vs = c(local_tree_vs, flatten_eutree(entry))
      }
    }
  }
  return(local_tree_vs)
}


#' Find match for a tree tip label in db
#'@export
match_tip_label = function(name, db, treedb){

  # match to hash
  if (stringr::str_sub(name, -1,-1) != 'h' & (stringr::str_length(rev(stringr::str_split(name, 'h')[[1]])[[1]]) == 8)  ) {

    hash = rev(stringr::str_split(name, 'h')[[1]])[[1]]
    match = match_ref(R = list(H = hash, N = 'ignore'), db, ignore_name = T)

    if (!is.null(match)) return(match)
  }


  # match by name
  base_name = stringr::str_split(name, '/')[[1]][2:4]
  base_name[[3]] = stringr::str_split(base_name[[3]], '_')[[1]][[1]]
  base_name_underscore = paste0(base_name, collapse =  '/')
  base_name = stringr::str_replace(base_name_underscore, '_', ' ')

  # get all viruses from db containing the root virus name
  db_names = sapply(db, function(v)v$N)
  db_idx = which(stringr::str_detect(db_names, base_name) | stringr::str_detect(db_names, base_name_underscore))
  print(name)
  print(db_idx)
  base_name_match_dbs = db[which(stringr::str_detect(db_names, base_name) | stringr::str_detect(db_names, base_name_underscore))  ]


  #' get entry from treedb
  treedb_names = sapply(treedb, function(v)v$n)
  base_name_match_treedb = treedb[[which(stringr::str_detect(treedb_names, name))  ]]



  # match to hi name from treedb
  if ('h' %in% names(base_name_match_treedb)){
    hi_names = base_name_match_treedb$h
    for (base_name_match_db in base_name_match_dbs){
      for (seq_entry in base_name_match_db$s){
        if (!is.null(seq_entry$h) & any(seq_entry$h %in% hi_names)){
          seq_entry_seqadd = with_sequence.seq_entry(seq_entry, list(base_name_match_treedb), db)
          if (stringr::str_detect(seq_entry_seqadd$a, base_name_match_treedb$a)){
            return(with_sequence.seq_entry(seq_entry, list(base_name_match_treedb), db))

          }
        }
      }
    }
  }

  # direct sequence match
  else{
    for (base_name_match_db in base_name_match_dbs){
      for (seq_entry in base_name_match_db$s) {
        seq_entry_seqadd = with_sequence.seq_entry(seq_entry, list(base_name_match_db), db)
        if (stringr::str_detect(seq_entry_seqadd$a, base_name_match_treedb$a)) return(seq_entry_seqadd)

      }
    }
  }

  warning('No match for: ', name)
  return(list(c = list('Clade missing')))

}
