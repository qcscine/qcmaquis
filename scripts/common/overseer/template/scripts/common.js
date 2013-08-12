var rank;

function latchLocks(){
    if(lockaction == true) return false;
    lockaction = true;
    return true;
}
function clear(){
    lockaction = false;
}
function findId(id, container){
    var match;
    container.children().each(function(){
	if($(this).data('id') == id) match = $(this);
    });
    return match;
}

function get_rank(){

    var filename = window.location.href.substr(window.location.href.lastIndexOf("/")+1);
    rank = filename.substr(filename.indexOf(".")+1);
    rank = rank.substr(0, rank.indexOf("."));

    if(!rank){
        if(window.location.href.lastIndexOf("#") > 1) 
            rank = window.location.href.substr(window.location.href.lastIndexOf("#")+1);
        else{
            rank = $("div.selector:first").text();
            rank = rank.substr(rank.indexOf(".")+1);
        }
        $("div.selector").each(function(){
            var i = $(this).text();
            i = i.substr(i.indexOf(".")+1);
            $(this).html("<a class='link' href=\"#"+i+"\">"+i+"</a>");
        }); 
    }else{


    }
}

function adjust(){
    var count = 0;
    $("div.region").each(function(){
        var title = $(this).data("title");
        $(this).html("<a href='log."+rank+"."+count+".html'>"+title+"</a>");
        count++;
    });
    var offset = 50;
    var left = ($(window).width() - ($("div.selector").size()+1) * offset)/2
    $("div.selector").each(function(){
        var link = $(this).find("a.link").each(function(){
            if(rank == $(this).text()) $(this).css({color: "red"});
            else $(this).css({color: "black"});
        });
        $(this).css({left: left});
        left += offset;
    });
}

$(document).ready(function() {

    get_rank(); adjust();

    $("a.link").click(function(){
        rank = $(this).text(); adjust();
    });

    $("div.rename").each(function(){
        var object_id = $(this).data("operand");
        var object_label = $(this).data("to");
        if($(this).data("to") == $(this).data("from")){ $(this).remove(); return; }
        $(this).html("<small>rename:</small> "+object_label+"["+object_id+"]");
    });

});
