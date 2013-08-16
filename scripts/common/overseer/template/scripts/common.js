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
        if($(this).find("a.index").length) $(this).find("a.index").attr("href", "log."+rank+"."+count+".html");
        else $(this).html("<a class='index' href='log."+rank+"."+count+".html' target='screen'>"+title+"</a>");
        if($(this).find("a.index").css("color") == "rgb(255, 0, 0)"){
            $("iframe#screen").attr("src", "log."+rank+"."+count+".html");
        }
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
    $("body").append("<iframe id=\"screen\"></iframe>");
    if($("div.selector").length){ 
        $("body").append("<div id='selector-container'></div>"); 
        $('div#selector-container').append($("div.selector"));
        $("body").prepend("<div id='region-container'></div>"); 
        $("div.region").appendTo($("div#region-container"));
        $("a.index").live("click", function(){
            $("a.index").css({color: "black"});
            $(this).css({color: "red"});
        });
    }

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

    $(window).keydown(function(e) {
        var keyCode = e.keyCode || e.which; 
        if (keyCode == 9){ 
            e.preventDefault();
            var current;
            $("div.selector").each(function(){
                if(rank == $(this).find("a.link").text()) current = $(this);
            });
            var neighbor = current.next();
            if(!neighbor.length) neighbor = $("div.selector:first");
            neighbor.find("a.link").click();
        } 
    });
});
